
#ifndef CEST_MODELS_H
#define CEST_MODELS_H

#ifdef _WIN32
#ifdef fabber_cest_EXPORTS
#define FABBER_CEST_API __declspec(dllexport)
#else
#define FABBER_CEST_API __declspec(dllimport)
#endif
#define CALL __stdcall
#else
#define FABBER_CEST_API
#define CALL
#endif

#include "fabber_core/fwdmodel.h"

extern "C"
{
    FABBER_CEST_API int CALL get_num_models();
    FABBER_CEST_API const char *CALL get_model_name(int index);
    FABBER_CEST_API NewInstanceFptr CALL get_new_instance_func(const char *name);
}

#endif
