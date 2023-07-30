#Cameron's code for job scheduling
#This code is in the script directory, the init.sh is one level above

import subprocess
import os.path
import time

#TODO: formatting the executables based on use inputs
def get_exe_syntax(logfile, scriptname: str):
    with open(logfile, 'r') as f:
        last_line = f.readlines()[-1]
        opts = last_line.split(' ')
        return [scriptname, 
                '-d', opts[1], 
                '-n', opts[3], 
                '-r', opts[5]]

#INI_EXE = get_exe_syntax('../.log.txt', '../input.sh')

QUEUE = [["EXECUTABLE", "PATH TO SCRIPT OUTPUT"],
         ["EXECUTABLE", None]]

SLEEP_TIME = 60*5

def main():
    for process in QUEUE:
        subprocess.run(process[0], shell=False)
        print('Ran script: ' + process[0])

        if process[1] is None: # Skip waiting for jobs without an output
            return

        while True:
            if os.path.isfile(process[1]):
                print('Found output file: ' + process[1])
                #TODO: check whether the file is finished being written to
                #TODO: input format validation
                
                break
            else:
                print('Waiting for output file: ' + process[1])
                time.sleep(SLEEP_TIME)

if __name__ == "__main__":
    main()

# EOF