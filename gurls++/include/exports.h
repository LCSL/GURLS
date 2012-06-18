#ifndef GURLS_EXP_H
#define GURLS_EXP_H

#ifdef _WIN32
#  ifdef _GURLS_STATIC
#    define GURLS_EXPORT __stdcall
#  elif defined _GURLS_EXPORTS
#    define GURLS_EXPORT  __declspec( dllexport )
#  else
#    define GURLS_EXPORT  __declspec( dllimport )
#  endif
#else
#  define GURLS_EXPORT
#endif

#endif // GURLS_EXP_H