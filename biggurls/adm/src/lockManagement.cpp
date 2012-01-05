#include "mex.h"
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/unistd.h>

#include "../include/retCode.h"
#include "../include/lockLib.h"

int acquireLock(int fLock)
{
		struct flock lock;
	
		lock.l_type 	= F_WRLCK;
		lock.l_whence 	= SEEK_SET;
		lock.l_start 	= 0;
		lock.l_len 		= 0;

		if(fcntl(fLock, F_SETLK, &lock) == -1)
		{
//				mexPrintf("Could not set lock\n");
				return LOCK_FAIL;
		}

		/*
		 *  At this point we have exclusive access to the 
		 *  state file, we are going to read what is left to do,
		 *  state that we are going to do it, clean up and leave.
		 */

		
		return LOCK_SUCCESS;
}

int releaseLock(int fLock)
{
		struct flock lock;
	
		lock.l_type 	= F_UNLCK;
		lock.l_whence 	= SEEK_SET;
		lock.l_start 	= 0;
		lock.l_len 		= 0;

		if(fcntl(fLock, F_SETLK, &lock) == -1)
		{
//				mexPrintf("Could not remove lock\n");
				return LOCK_FAIL;
		}
		return LOCK_SUCCESS;
}


