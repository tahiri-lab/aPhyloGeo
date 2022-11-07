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
        self.resultList = Manager().list()
        self.function = function

        self.processlist= []
        self.mem1 = 1  
        self.memA = psutil.virtual_memory()[1]
        self.nbAllowed = 1
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
        self.terminalUpdate()
        self.tasks.value+= 1
        self.started.value +=1

        self.resultList.append(self.function(arg))
        self.mem1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000
        self.tasks.value -= 1
        self.finished.value += 1
        self.terminalUpdate()


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
        print("Starting calibration run, this might some minutes")
        self.startTime = time.time()

        #self.executeOnce(self.args.pop(0))  # Calibration run; needs to be alone to set the required memory usage of a unique task
   
        p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs
        p.start()
        time.sleep(1)
        while(True):
            if(self.tasks.value==0):
                self.memUpdate()
                p.join()
                break
            else:
                current_process = psutil.Process(os.getpid())
                mem = current_process.memory_percent()
                for child in current_process.children(recursive=True):
                    mem+= child.memory_percent()
                mem =(self.memA/100*mem)/1000000000
                if (self.mem1<mem):
                    self.mem1=mem

                #print("*****\n",mem)
                self.terminalUpdate()
                
                time.sleep(0.1)
        
        while(len(self.args)!=0):
            if self.tasks.value < self.nbAllowed:
                p = Process(target=self.executeOnce, args=([self.args.pop(0)])) #Multiprocess runs
                
                self.processlist.append(p)
                p.start()
                time.sleep(0.1)
            else:
                time.sleep(0.1)
            self.terminalUpdate()

        while(True):
            if(self.tasks.value==0):
                for p in self.processlist:
                    p.join()
                break
            else:
                self.terminalUpdate()

        return self.resultList

    """
    Method that sets the baseline for memory calculation
    + output some information to the terminal
    All memory values are in bytes

    Variables:
        memBuffer  (double) Amount of byte to substract from the available RAM for safety purposes
    """
    def memUpdate(self):
        memBuffer = 1000000000
        #self.mem1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000
        self.nbAllowed = math.floor((((self.memA-memBuffer)/self.mem1)/1000000000))

        print("Available memory: ", round(self.memA/1000000)/1000, "Gb"+
                    "                                                         ")
        print("needed memory:    ",round(self.mem1*1000)/1000, "Gb"+
                    "                                                         ")
        print("Max processes:    ", self.nbAllowed, end="\n")
    
    """
    Method that constantly updates the user about the currently run tasks
    """
    def terminalUpdate(self):
        s=self.started.value
        a=self.amount
        f=self.finished.value
        nowTime = time.time()
        eTime = round((nowTime - self.startTime)*100)/100
        print(" {}/{} Started    {}/{} Finished    Time elapsed: {} seconds".format(s,a,f,s,eTime), end="\r")
    
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
        self.terminalUpdate()
    
        for a in range(len(self.args)):
            p = Process(target=self.executeSmall, args=([self.args.pop(0)])) #Multiprocess runs
            self.processlist.append(p)
            p.start()
            self.terminalUpdate()

        for p in self.processlist:
            p.join()
            self.terminalUpdate()

        self.terminalUpdate()

        return self.resultList


