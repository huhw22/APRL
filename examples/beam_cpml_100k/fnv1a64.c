/* Generate the native loader's FNV-1a identifiers when assembling the dataset. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    unsigned char buffer[65536];
    for (int i = 1; i < argc; ++i) {
        FILE *in = fopen(argv[i], "rb");
        if (!in) { perror(argv[i]); return 1; }
        uint64_t state = UINT64_C(14695981039346656037);
        size_t n;
        while ((n = fread(buffer, 1, sizeof buffer, in)) > 0) {
            for (size_t j = 0; j < n; ++j)
                state = (state ^ buffer[j]) * UINT64_C(1099511628211);
        }
        if (ferror(in)) { perror(argv[i]); fclose(in); return 1; }
        fclose(in);
        printf("%016" PRIx64 "\n", state);
    }
    return 0;
}
