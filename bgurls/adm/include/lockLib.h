#ifndef _LOCKLIB_H_
#define _LOCKLIB_H_


int acquireLock(int fLock);
int releaseLock(int fLock);

int getWork(const char *stateFileName, int fLock);
int reportWork(const char *stateFileName, int fLock, int finishedBlockID);



#endif
