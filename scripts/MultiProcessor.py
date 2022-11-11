from multiprocess import Process,Value,Manager
import psutil
import time
import resource
import math
import os

"""
Class that harnesses the power of multiprocessing and package for ease of use. as a single line.
It is mainly aimed at very many large tasks that could be run on supercumputers with insane amount of RAM.

It's methods take in a list of arguments and a method. 
Each position on the first level of the list will be the arguments given to the child process running the given method

The methods used must only take a list as input; and then, deal with it<s division into other variables

ex.:
    def g(args):
        n1=args[0]
        n2=args[1]
        return n1+n2
    
    list = [[0,1],[1,1],[2,1]]

    mps = Multi(list,g).processingSmallData()

    >>[0,1,2]
"""
class Multi:

    """
    Constructor of the multiprocessing tool 

    Args:
        args        (list)      The list of value to run. One process will be born per 1st level member
        function    (method)    The method that needs to run as multiple instance

    Variables:
        args        (Manager().list()) A list accessible from child processes; for the arguments
        resultList  (Manager().list()) A list accessible from child processes; for the return values
        function    (method)    The method that needs to run as multiple instance

        processlist (list)      All of the child processes    
        mem1        (double)    The memory needed to run a single child
        memA        (double)    The current available memory
        nbAllowed   (int)       The amount of theoretical child processes that can fit in memory
        tasks       ()          The amount of child processes currently executing
        original    (String)    The DNA sequence to compare to

        started     (Value(i))  Number of started processes
        finished    (Value(i))  Number of finished processes
        amout       (int)       Total of all processes
    """
    def __init__(self,args,function):


        self.args = Manager().list(args)
        self.processes = Manager().list()
        self.resultList = Manager().list()
        self.function = function

        self.processlist= []
        self.mem1 = Value("f",1)  
        self.memA = Value("f",1)
        self.memM = Value("f",psutil.virtual_memory()[1])
        self.nbAllowed = Value('i',1)
        self.tasks = Value('i',0)

        self.started = Value('i',0)
        self.finished = Value('i',0)
        self.amount = len(args)

        self.startTime = 0
    
    """
    The method executed by the RAM hungry processes

    Return:
        Nothing, but the return value of the executed method is passed to a global multiprocessing-friendly list
    """
    def executeOnce(self,arg):
        self.tasks.value+= 1
        self.started.value +=1
        self.processes.append(os.getpid())

        self.resultList.append(self.function(arg))
        self.mem1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000
        
        self.tasks.value -= 1
        self.finished.value += 1


    """
    Method for executing mutliprocessing on tasks that demand a large amount of individual memory.
    Will start as many child processes as the available RAM permits

    Variables:
        p   (Process)   Representes a single child process

    Return:
        The multiprocess-friendly list, that was updated by each child

    Errors:
        If other application reduce the avalable RAM mid-execution,
        Multiprocess outputs "Killed" and kills the child.

    """
    def processingLargeData(self):
        print("Starting calibration run, the first instance runs alone")
        self.startTime = time.time()

        t = Process(target=self.terminalUpdate)
        t.start()

        p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs
        p.start()

        self.processes.append(os.getpid())
        m = Process(target=self.memUpdate)
        m.start()
        
        time.sleep(1)
        while(self.tasks.value>0):
            #WAIT FOR CALLIBRATION
            time.sleep(0.1)
        
        while(len(self.args)!=0):
            if self.tasks.value < self.nbAllowed.value:
                p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs
                
                self.processlist.append(p)
                p.start()
                time.sleep(0.1)
            else:
                time.sleep(0.1)


        for p in self.processlist:
            p.join()
        time.sleep(1)
            
        t.terminate()
        m.terminate()
        print("\033[F", end="")
        print("\033[B"*5, flush=True)
        print("completed with ",str(self.amount-self.finished.value)," errors")

        return self.resultList

    """
    Method that sets the baseline for memory calculation
    + output some information to the terminal
    All memory values are in bytes

    Variables:
        memBuffer  (double) Amount of byte to substract from the available RAM for safety purposes
    """
    def memUpdate(self):
        memBuffer = 5 * 1000000000
        while True:
            self.memA.value = psutil.virtual_memory()[1]
            for child in self.processes:
                try:
                    mem = psutil.Process(child).memory_percent()*self.memM.value/100000000000
                    if (self.mem1.value<mem):
                        self.mem1.value=mem
                except:
                    ""

            self.nbAllowed.value = math.floor((((self.memM.value-memBuffer)/self.mem1.value)/1000000000))
            if self.nbAllowed.value<1:
                self.nbAllowed.value=1
    
    """
    Method that constantly updates the user about the currently run tasks
    """
    def terminalUpdate(self):
        while True:
            s=self.started.value
            a=self.amount
            f=self.finished.value
            nowTime = time.time()
            eTime = round((nowTime - self.startTime)*100)/100
            print("Maximum memory:   ", round(self.memM.value/1000000)/1000, "Gb        ", end="\n", flush=True)
            print("Available memory: ", round(self.memA.value/1000000)/1000, "Gb        ", end="\n", flush=True)
            print("needed memory:    ",round(self.mem1.value*1000)/1000, "Gb        ", end="\n", flush=True)
            print("Active processes: ", str(s-f) + "/" + str(self.nbAllowed.value) + "            ", end="\n", flush=True)
            print(" {}/{} Started    {}/{} Finished    Time elapsed: {} seconds".format(s,a,f,s,eTime), flush=True)
            print("\033[F", end="",flush=True)
            print("\033[A"*5, flush=True)
            time.sleep(0.01)
            
    
    """
    The method executed by small processes

    Return:
        Nothing, but the return value of the executed method is passed to a global multiprocessing-friendly list
    """
    def executeSmall(self,arg):
        self.started.value +=1
        self.startTime = time.time()

        result = self.function(arg)
        self.resultList.append(result)

        self.finished.value += 1
    
    """
    Method for executing mutliprocessing on tasks that demand little to no memory.
    Will immedialty start all the child processes
    The packaging causes some marginal time lost; 
    Only use for methods that take at least a second to run
        below that, a for loop is likely much faster

    Variables:
        p   (Process)   Representes a single child process
        a   (None)      Exists only to permet the for loop

    Return:
        The multiprocess-friendly list, that was updated by each child

    Errors:
        If other application reduce the avalable RAM mid-execution,
        Multiprocess outputs "Killed" and kills the child.

    """
    def processingSmallData(self):
        self.startTime = time.time()
        t = Process(target=self.terminalUpdate).start()

    
        for a in range(len(self.args)):
            p = Process(target=self.executeSmall, args=([self.args.pop(0)])) #Multiprocess runs
            self.processlist.append(p)
            p.start()

        for p in self.processlist:
            p.join()

        t.terminate()

        return self.resultList


