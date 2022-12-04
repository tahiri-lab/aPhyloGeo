from multiprocess import Process,Value,Manager
import psutil
import time
import resource
import math
import os


class Multi:
    """
    Class that harnesses the power of "multiprocess" and wraps it, for ease of use, as a single callable line.
    It is mainly aimed at very many large tasks that could be run on supercomputers with an ridiculous amount of RAM.
    """

    def __init__(self,args,function):
        """
        Constructor of the multiprocessing tool 

        Args:
            args        (list)      The list of value to run. One process will be born per 1st level member
            function    (method)    The method that needs to run as multiple instance

        The constructor needs two thing, a method, and a list of lists. 
        Each sub-list of the list will be the list of arguments given to the child process running the given method
        ***The methods used must only take a single list as input; and then, deal with it's division into other variables***

        ex.: Run 3 processes of g()

            def g(args):
                n1=args[0]
                n2=args[1]
                return n1+n2
            
            list = [ [0,1], [1,1] ,[2,1] ]

            mps = Multi( list, g ).processingSmallData()        #note that "g" is not written "g()" we want to pass it, not run it
            print( mps )

            >>[1,2,3]

        Variables:
            args        Manager().list()    A list accessible from child processes; for the arguments
            processes   Manager().list()    A list accessible from child processes; list of all the process ID started in this object
            resultList  Manager().list()    A list accessible from child processes; for the return values
            function    method              The method that needs to run as multiple instances

            processlist list        All of the child processes    
            mem1        Value(f)    A value accessible from child processes; The memory needed to run a single child
            memA        Value(f)    A value accessible from child processes; The current available memory
            memT        Value(f)    A value accessible from child processes; The total amount of available memory at the creation of the object
            nbAllowed   int         The amount of theoretical child processes that can fit in FREE memory
            maxAllowed  int         The amount of theoretical child processes that can fit in TOTAL memory
            tasks       Value(i)    A value accessible from child processes; The amount of child processes currently running

            started     Value(i)    A value accessible from child processes; Number of started processes
            finished    Value(i)    A value accessible from child processes; Number of finished processes
            amount      int         The amount of processes that the object will start

            startTime   float       the present time
            timeForOne  Value(f)    A value accessible from child processes; Time it takes for the first iteration to complete
            rewrite     dict        ints used the the terminal update rewriting process; amount of lines rewinded if large or small
        """

        self.args = Manager().list(args)
        self.processes = Manager().list()
        self.resultList = Manager().list()
        self.function = function

        self.processlist= []
        self.mem1 = Value("f",1)  
        self.memA = Value("f",1)
        self.memT = Value("f",psutil.virtual_memory()[1])
        self.nbAllowed = Value('i',1)
        self.maxAllowed = Value('i',1)
        self.tasks = Value('i',0)

        self.started = Value('i',0)
        self.finished = Value('i',0)
        self.amount = len(args)

        self.startTime = 0
        self.timeForOne = Value("f",0)
        self.rewrite= {True:11,False:6}
 
    def executeOnce(self,arg):
        """
        The method that is ran as a single process

        Args:
            args    list    The list of arguments given to the method

        Return:
            Nothing None    But the return value of the executed method is passed to self.resultList
        """
        self.tasks.value+= 1
        self.started.value +=1
        self.processes.append(os.getpid())

        self.resultList.append(self.function(arg)) #execution of the passed method
        
        self.tasks.value -= 1
        self.finished.value += 1

    def processingLargeData(self):
        """
        Method for executing mutliprocess on tasks that demand a LARGE amount of individual memory.
        Will, first, run a single process then, will start as many child processes as the available RAM permits
        Starting new ones as the RAM is freed

        Variables:
            p   (Process)   Representes a single child process

        Return:
            The multiprocess-friendly list, that was updated by each child

        Errors:
            If other application reduce the avalable RAM mid-execution,
            Multiprocess outputs "Killed" and kills the child.

        """
        print("    Starting multiprocessing, this might take some time\n    The first process is ran alone for calibration purposes")
        self.startTime = time.time()

        p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs once alone
        p.start()

        print("\033[B"*self.rewrite[False], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        self.processes.append(os.getpid())              #adds the main thread on the list of processes to keep an eye on

        alfred = Process(target=self.buttler,args=([True]))  #ask the buttler to start complimentary processes
        alfred.start()
        
        time.sleep(1)   #give it a second to open the process, so it doesn't skip the while()
        while(self.tasks.value>0):
            #wait for the calibration process to finish
            time.sleep(0.1)
        self.timeForOne.value=time.time()-self.startTime
        
        while(len(self.args)!=0):
            if (self.tasks.value < self.maxAllowed.value) & (self.nbAllowed.value >= 1):
                p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs the rest of the processes
                
                self.processlist.append(p)
                p.start()
                time.sleep(0.1)
            else:
                time.sleep(0.1)


        for p in self.processlist:
            p.join()
        time.sleep(1)   #give it a second to close the processes, do not remove

        alfred.terminate()   #sorry Alfred, we have to let you go
        #These weird prints need to be done because there is no telling where the terminal is at termination time, better make sure it's clean
        print("\r", end="")
        print("\033[B"*self.rewrite[True], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        print("Completed with ",str(self.amount-self.finished.value)," errors\n") #if a process was killed or didn't finished; it will be know here

        return self.resultList
  
    def buttler(self,memBloc):
        """
        Ran as a child process, the buttler will constantly run other methods forever.
        In this case, it:
            updates de memory capacity and 
            prints updates on the terminal.
        It exists so not to bottleneck the main thread.

        Uses timers to execute it's methods because time.sleep() it processor hungry if constantly called
        """
        terminal = time.time()
        mem = time.time()
        while True:
            now = time.time()
            if now-terminal>0.1:    #has 0.1 second passed?
                self.terminalUpdate(memBloc)
                terminal = now
            if now-mem>1 & memBloc: #has 1 second passed?
                self.memUpdate()
                mem = now

    def memUpdate(self):
        """
        Method that sets the baseline for memory calculation
        + output some information to the terminal
        All memory values are in bytes

        This method is ran from the buttler() and updates every second

        Variables:
            memBuffer  double   %Amount of bytes to substract from the available RAM for safety purposes
            mem        double   Amount of bytes 
        """

        memBuffer = 0.9 #%
        self.memA.value = psutil.virtual_memory()[1]*memBuffer
        for child in self.processes:
            try:    #in a try/except because processes ID are never removed from the list
                mem = psutil.Process(child).memory_full_info()[8] #uss memory usage; humch much is this process using NOW
                if (self.mem1.value<mem):
                    self.mem1.value=mem #does it for the whole run in case this maximum is increased by future childs
            except:
                ""

        self.nbAllowed.value = math.floor(((self.memA.value/self.mem1.value)))
        if self.nbAllowed.value < 1:
            self.nbAllowed.value = 1    #Need to at least be able to start a single process
        self.maxAllowed.value = math.floor(((self.memT.value/self.mem1.value)))
        if self.maxAllowed.value < 1:
            self.maxAllowed.value = 1   #Need to at least be able to start a single process

    def terminalUpdate(self, memBlock):
        """
        Method that constantly updates the user about the currently run tasks

        This method is ran from the buttler() and updates every 0.1 second
        """

        s=self.started.value
        a=self.amount
        f=self.finished.value
        nowTime = time.time()
        eTime = round((nowTime - self.startTime)*10)/10 #current execution time
        
        print("---")
        if(memBlock): #block of prints used only by processingLargeData()
            print("Available memory: ", round(self.memA.value/10000000)/100,"/",round(self.memT.value/10000000)/100, "Gb        ", end="\n", flush=True)
            print("Active processes: ", str(self.tasks.value) + " / " + str(self.maxAllowed.value) + "            ", end="\n", flush=True)
            print("Min memory per:   ",round(self.mem1.value/10000000)/100, "Gb        ", end="\n", flush=True)
            print("Time for one:     " ,round(self.timeForOne.value*10)/10," seconds               ", flush=True)
            print("---")
        print("Started:          ",s,"/",a,"   ",round(s/a*100),"%           ", flush=True)
        print("Finished:         ",f,"/",s,"   ",round(f/a*100),"%           ", flush=True)
        print("Time elapsed:     " ,eTime," seconds               ", flush=True)
        print("---")

        print("\r", end="",flush=True)
        print("\033[A"*self.rewrite[memBlock], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up

    def executeSmall(self,arg):
        """
        The method executed by processingSmallData
        
        Return:
            Nothing, but the return value of the executed method is passed to a global multiprocessing-friendly list
        """
        self.started.value +=1

        result = self.function(arg)
        self.resultList.append(result)

        self.finished.value += 1
    
    def processingSmallData(self):
        """
        Method for executing mutliprocess on tasks that demand little to no memory.
        Will immedialty start all the child processes
        The packaging causes some marginal time lost; 
        Only use for methods that take at least a second to run
            below that, a for loop is likely much faster

        Variables:
            p   Process   Representes a single child process
            a   None      Exists only to permit the for loop

        Return:
            The multiprocess-friendly list, that was updated by each child

        Errors:
            If other application reduce the avalable RAM mid-execution,
            Multiprocess outputs "Killed" and kills the child.

        """

        self.startTime = time.time()
        alfred = Process(target=self.buttler,args=([False]))
        alfred.start()

    
        for a in range(len(self.args)):
            p = Process(target=self.executeSmall, args=([self.args.pop(0)])) #Multiprocess runs
            self.processlist.append(p)
            p.start()

        for p in self.processlist:
            p.join()

        time.sleep(1) #wait for the processes to close; do not remove
        finishedTime = round(time.time()-self.startTime)*10/10

        alfred.terminate()
        print("\r", end="")
        print("\033[B"*self.rewrite[False], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        print("Completed ",len(self.resultList),"tasks in ",finishedTime,"seconds")
 
        return self.resultList
