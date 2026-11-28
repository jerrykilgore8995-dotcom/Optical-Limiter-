#pragma once

// Specialize each generated stamp for its fixed terminals and model constants.
#if defined(_MSC_VER) && !defined(__clang__)
#define OPTICAL_STAMP_INLINE inline
#elif defined(_MSC_VER)
#define OPTICAL_STAMP_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define OPTICAL_STAMP_INLINE inline __attribute__((always_inline))
#else
#define OPTICAL_STAMP_INLINE inline
#endif
