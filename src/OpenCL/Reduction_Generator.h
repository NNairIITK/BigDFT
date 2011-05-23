#ifndef REDUCTION_GENERATOR_H
#define REDUCTION_GENERATOR_H

#ifdef __cplusplus
extern "C" const char* generate_reduction_program(struct bigdft_device_infos * infos);
#else
const char* generate_reduction_program(struct bigdft_device_infos * infos);
#endif

#endif
