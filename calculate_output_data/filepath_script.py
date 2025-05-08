import os


# Максимальное значение для int32
# INT32_MAX = 2**31 - 1
# INT32_MIN = -2**31 - 1

# print(INT32_MAX)

# 2100000000
# 2150000000

# num = -2144967296

# result = num % (2**32)

num_max = 2**32

# result = INT32_MAX + (-2144967296 - INT32_MIN)

# print(result)

# Путь к директории
directory_path = 'calculate_output_data/Lx20.0_Ly20.0_Nx64_Ny64_eps0.40_dt1.0e-08_Periodic_2d_half_from_file_-0.9'

numbers = []
# Проход по всем файлам в директории
for filename in os.listdir(directory_path):
    if filename in ('params.yaml', 'energies.bin'):
        continue
    file_number = int(filename.removesuffix('.bin'))
    
    if file_number % 10 != 0:
        file_number = file_number % num_max

        if file_number % 10 != 0:
            file_number += num_max

    old_file_path = os.path.join(directory_path, filename)
    new_file_path = os.path.join(directory_path, f'{file_number}.bin')

    os.rename(old_file_path, new_file_path)
    print(f'Renamed: {old_file_path} -> {new_file_path}')


