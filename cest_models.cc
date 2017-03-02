/* dualecho_models.cc Shared library functions for dualecho models

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "fwdmodel_cest.h"

extern "C" {
int get_num_models()
{
    return 1;
}

const char *get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "cest";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr get_new_instance_func(const char *name)
{
    if (string(name) == "cest")
    {
        return CESTFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
