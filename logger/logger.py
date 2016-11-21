#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
    
from time import time


# Print iterations progress
def printProgress (iteration, total, prefix = 'Progress: ', suffix = '', decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    fmt             = "{0:." + str(decimals) + "f}"
    percents        = fmt.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r{} |{}| {}%% {}'.format(prefix, bar, percents, suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

class ProgressBar:
    def __init__(self, prefix = 'Progress:', suffix = '', decimals = 1, barLength = 50, nUpdates=50, start=False):
        self.t0 = time()
        self.prefix = prefix
        self.suffix = suffix
        self.barLength = barLength
        self.fmt = "{0:." + str(decimals) + "f}"
        self.nUpdates = nUpdates
        self.update = -100*sys.float_info.epsilon
        if start:
            self.show(0)
        
    def show(self,fraction):
        if fraction > self.update:
            self.update += 1/self.nUpdates
            percents        = self.fmt.format(100 * fraction)
            filledLength    = int(round( self.barLength * fraction))
            bar             = '█' * filledLength + '-' * (self.barLength - filledLength)
            sys.stdout.write('\r{} |{}| {}%% {}'.format(self.prefix, bar, percents, self.suffix))
            if fraction >= 1:
                sys.stdout.write(' in {:.3g} s.'.format(time()-self.t0))
                self.update = 0 # flag that we are done
            else:
                while fraction>self.update:
                    self.update += 1/self.nUpdates
            sys.stdout.flush()         
        
    def finish(self):
        if self.update!=0:
            self.show(1)
    
    def __del__(self):
        self.finish()
        
    def output(self,*arg):
        # print cannot be defined in python 2.7
        s = '\n'
        for a in arg:
            s+=str(a)+' '
        s += '\n'
        sys.stdout.write(s)
        sys.stdout.flush()

class Log():
    """
    context manager class for logging
    """
    level = 0
#     sth_printed_before_exit = False
    current_logger=None
    out = None
    
    @staticmethod
    def Reset():
        Log.level = 0
#         Log.sth_printed_before_exit = False
        Log.current_logger = None
        
    @staticmethod
    def Indent(msg):
        s='\n'
        if Log.level==0:
            return s+msg
        else:
            return s + ('|'.ljust(3))*Log.level + msg
    
    @staticmethod
    def Output(*arg):
        if Log.current_logger:
            Log.current_logger.output(*arg)
        else:
            # equivalent to print(*arg) but works in python 2.7 too
            for a in arg:
                sys.stdout.write(str(a)+' ')            
            sys.stdout.write('\n')
            sys.stdout.flush()
            
    @staticmethod
    def Show_progress(fraction):
        Log.current_logger.show_progress(fraction)
    
    @staticmethod
    def Print_final(*arg):
        Log.current_logger.output_final(*arg)
    
    def __init__(self, title, progressBar=None, out=None):
        """
        """
        if Log.level==0:     
            if out:
                if isinstance(out,str):
                    Log.out = open(out,'w')
                else: #hopefully a file object (but this is not of type file in python 3 
                    Log.out = out
        self.sth_printed_before_exit = (not progressBar is None)
        if Log.current_logger:
            Log.current_logger.sth_printed_before_exit = True # previous level logger, not this onw
            
        self.title = title
        self.progressBar = progressBar
        s = Log.Indent(title)+' ... '
        sys.stdout.write(s)
        if Log.out:
            Log.out.write(s)
        Log.level += 1
        self.level = Log.level
        if self.progressBar:
            self.progressBar.prefix = Log.Indent(self.progressBar.prefix)[1:]
            sys.stdout.write('\n')
        self.final_message = ''
        
        Log.current_logger = self
        
    def __enter__(self):
        self.t0 = time()
        if self.progressBar:
            self.progressBar.t0 = self.t0
        return self
    
    def __exit__(self, *args):
        if self.progressBar:
            self.progressBar.finish()
            
        if self.final_message:
            sys.stdout.write(self.final_message)
            self.sth_printed_before_exit = True
            if Log.out:
                Log.out.write(self.final_message)
        
        Log.level-=1
        
        
        s = 'completed in {:.3g} s.'.format(time()-self.t0)
        if self.sth_printed_before_exit:
            s = Log.Indent(self.title+': '+s)
        sys.stdout.write(s)        
        if Log.out:
            Log.out.write(s)
        
        if Log.level==0:
            Log.Reset()
            self.output('\n')

    def show_progress(self,fraction):
        if self.progressBar:
            if self.level != Log.level:
                #silent removal of self.ProgressBar
                self.progressBar.update=0
                self.progressBar = None
            else:
                self.progressBar.show(fraction)

    def indent(self,*arg):
        s = Log.Indent('')
        if arg:
            s += str(arg[0])
            for a in arg[1:]:
                s += ' '+str(a)
        return s
                    
    def output(self,*arg):
        self.sth_printed_before_exit = True
        s = self.indent(*arg)
        sys.stdout.write(s)
        if Log.out:
            Log.out.write(s)
        if self.progressBar:
            sys.stdout.write('\n')
    
    def output_final(self,*arg):
        self.final_message += self.indent(*arg)
                    
#==================================================================================================
# test code below
#==================================================================================================
if __name__=='__main__':
    with Log('count to 1000') as logger:
        sum = 0
        for i in range(1000):
            sum +=i

    with Log('count to 1000') as logger:
        sum = 0
        for i in range(1000):
            sum +=i
        logger.output_final('hi, everything went fine!')
    
    with Log('count to 1000'):
        sum = 0
        for i in range(1000):
            if i==500:
                Log.Output('something happened at 500')
                Log.Output('yep')
            sum +=i
    
    with Log('nested count'):
        for j in range(5):
            with Log('count to 1000') as logger:
                sum = 0
                for i in range(1000):
                    if j==2 and i==500:
                        logger.output('something happened at 2:500')
                    sum +=i

    with Log('count to 1000',progressBar=ProgressBar()) as logger:
        sum = 0
        for i in range(1000):
            logger.show_progress((i+1)/1000)
            sum +=i

    with Log('nested count **'):
        for j in range(5):
            with Log('count to 1000',progressBar=ProgressBar()) as logger:
                sum = 0
                for i in range(1000):
                    Log.Show_progress((i+1)/1000)
                    if j==2 and i==500:
                        logger.output('something happened at',2,':',500)
                    sum +=i
                logger.output_final('hello, something happened ...')

    with Log('nested count with file',out='nested_count_with_file.txt'):
        for j in range(5):
            with Log('count to 1000',progressBar=ProgressBar()) as logger:
                sum = 0
                for i in range(1000):
                    Log.Show_progress((i+1)/1000)
                    if j==2 and i==500:
                        logger.output('something happened at',2,':',500)
                        logger.output_final('hello, something happened ...')
                    sum +=i

    from time import sleep
    
    # make a list
    items = list(range(0, 100))
    i     = 0
    l     = len(items)
    
    # Initial call to print 0% progress
    printProgress(i, l, prefix = 'printProgress:', suffix = 'Complete', barLength = 50)
    for item in items:
        # Do stuff...
        sleep(0.01)
        # Update Progress Bar
        i += 1
        printProgress(i, l, prefix = 'printProgress:', suffix = 'Complete', barLength = 50)

    p = ProgressBar(prefix='Just ProgressBar:')
    i=0
    for item in items:
        # Do stuff...
        sleep(0.01)
        # Update Progress Bar
        i += 1
        if i%20==0:
            p.output('something happened at',i)
        p.show(i/100)
    del p
    
    print('\n\ndone')
