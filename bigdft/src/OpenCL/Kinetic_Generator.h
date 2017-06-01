#ifndef KINETIC_GENERATOR_H
#define KINETIC_GENERATOR_H

#ifdef __cplusplus
extern "C" char* generate_kinetic_program(struct bigdft_device_infos * infos);
#else
char* generate_kinetic_program(struct bigdft_device_infos * infos);
#endif

#endif
