#ifndef BENCHMARK_GENERATOR_H
#define BENCHMARK_GENERATOR_H

#ifdef __cplusplus
extern "C" char* generate_benchmark_program(struct bigdft_device_infos * infos);
#else
char* generate_benchmark_program(struct bigdft_device_infos * infos);
#endif

#endif
