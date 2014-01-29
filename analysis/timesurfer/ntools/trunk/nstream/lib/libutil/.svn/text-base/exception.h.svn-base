#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

#include <sys/types.h>
#include <stdio.h>
#include <stdbool.h>

/* hack hack slash slash */

typedef struct Exception Exception;

struct Exception {
    int exception_number;
    char exception_text[512];
    char exception_format[512];
    bool enable_print;
};

#ifdef __cplusplus
extern "C"
{
#endif

Exception *Exception_create(const char *fmt, ...);
void Exception_destroy(Exception *ex);
void Exception_set(Exception *ex, const char *fmt, ...);
char *Exception_get(Exception *ex, const char *msg, int len);
void Exception_propagate(Exception *src, Exception *dst);
void Exception_propagate_with_text(Exception *src, Exception *dst, const char *text);
void Exception_enable_print(Exception *ex, bool print);
void Exception_printf(const char *fmt, ...);
void Exception_set_strerror(Exception *ex, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif

