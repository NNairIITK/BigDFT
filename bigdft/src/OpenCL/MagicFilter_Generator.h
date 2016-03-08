#ifndef MAGICFILTER_GENERATOR_H
#define MAGICFILTER_GENERATOR_H

#ifdef __cplusplus
extern "C" char* generate_magicfilter_program(struct bigdft_device_infos * infos);
#else
char* generate_magicfilter_program(struct bigdft_device_infos * infos);
#endif

#endif
