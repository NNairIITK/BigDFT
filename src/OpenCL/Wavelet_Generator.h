#ifndef WAVELET_GENERATOR_H
#define WAVELET_GENERATOR_H

#ifdef __cplusplus
extern "C" char* generate_ana_program(struct bigdft_device_infos * infos);
extern "C" char* generate_syn_program(struct bigdft_device_infos * infos);
#else
char* generate_ana_program(struct bigdft_device_infos * infos);
char* generate_syn_program(struct bigdft_device_infos * infos);
#endif

#endif
