#if !defined(_ALIGNMENT)
#define _ALIGNMENT 4096
#endif
#if !defined(_PADDING)
#define _PADDING 1
#endif

#if defined(_EXIT_ON_ERROR)
#define EXIT_ON_ERROR STOP
#else
#define EXIT_ON_ERROR RETURN
#endif
