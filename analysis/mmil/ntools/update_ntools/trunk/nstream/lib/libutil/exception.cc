#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include "exception.h"

Exception *Exception_create(const char *fmt, ...)
{
    Exception *ex;

    if ((ex = (Exception*)malloc(sizeof(Exception))) == NULL)
        return false;

    memset(ex, 0, sizeof(Exception));
    ex->enable_print = true;

    va_list args;
    va_start(args, fmt);
    vsnprintf(ex->exception_format, sizeof(ex->exception_format), fmt, args);
    va_end(args);

    return ex;
}

void Exception_destroy(Exception *ex)
{
    free(ex);
}

char *Exception_get(Exception *ex, const char *msg, int len)
{
    return strncpy((char*)msg, ex->exception_text, len);
}

void Exception_enable_print(Exception *ex, bool print)
{
    ex->enable_print = print;
} 

void Exception_set(Exception *ex, const char *fmt, ...)
{
    sprintf(ex->exception_text, "%s", ex->exception_format);

    va_list args;
    va_start(args, fmt);
    vsnprintf(ex->exception_text + strlen(ex->exception_text), sizeof(ex->exception_text) - strlen(ex->exception_text), fmt, args);
    va_end(args);
    if (ex->enable_print)
        Exception_printf("%s", ex->exception_text);
}

void Exception_set_strerror(Exception *ex, const char *fmt, ...)
{
    sprintf(ex->exception_text, "%s", ex->exception_format);

    va_list args;
    va_start(args, fmt);
    vsnprintf(ex->exception_text + strlen(ex->exception_text), sizeof(ex->exception_text) - strlen(ex->exception_text), fmt, args);
    va_end(args);

    snprintf(ex->exception_text + strlen(ex->exception_text), sizeof(ex->exception_text) - strlen(ex->exception_text), ": %s", strerror(errno));
    if (ex->enable_print)
        Exception_printf("%s", ex->exception_text);
}

void Exception_propagate(Exception *src, Exception *dst)
{
    Exception_set(dst, src->exception_text);
} 

void Exception_propagate_with_text(Exception *src, Exception *dst, const char *text)
{
    Exception_set(dst, "%s%s", text, src->exception_text);
}

void Exception_printf(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
//    if (EXCEPTION_OUTPUT_PIPE_fp != NULL)
//        vfprintf(EXCEPTION_OUTPUT_PIPE_fp, fmt, args);
//    else
        vfprintf(stderr, fmt, args);
    va_end(args);

//    if (EXCEPTION_OUTPUT_PIPE_fp != NULL)
//        fflush(EXCEPTION_OUTPUT_PIPE_fp);
}
    
