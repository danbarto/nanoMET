''' File lock mechanism for Cache
'''

# Standard imports
import os, errno, time

# Logger
import logging
logger = logging.getLogger(__name__)

# Wait until multithread process could acquire lock
def waitForLock(filename):
    lockAcquired = False
    while not lockAcquired:
      try:
           f = os.open(filename + "_lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY)
           os.close(f)
           lockAcquired = True
      except OSError as e:
           if e.errno == errno.EEXIST:  # Failed as the file already exists.
             time.sleep(1)
           else:  # Something unexpected went wrong
             logger.error( "Problem acquiring the lock" )
             exit(1)

def removeLock(filename):
    os.system("rm " + filename + "_lock")

