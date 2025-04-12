import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np


T_steps = None 
moment_in_frame = None
eps = None
N = None

with open("bin_data/params.txt") as file:
    T_steps, moment_in_frame = [int(file.readline()) for _ in range(2)]
    eps = float(file.readline()) 

frames = (T_steps // moment_in_frame) + 1 # Общее количество кадров
filenames = []


with open(f'bin_data/energies.bin', 'rb') as file:
    N = np.fromfile(file, dtype=np.int32, count=1)[0]
    energies = np.fromfile(file, dtype=np.float64, count=N)

print(energies)
x_extrems = []
y_extrems = []


phi_balance = []
for frame in range(frames):
    with open(f'bin_data/data_{frame * moment_in_frame}.bin', 'rb') as file:
        N = np.fromfile(file, dtype=np.int32, count=1)[0]
        phi = np.fromfile(file, dtype=np.float64, count=N)
    phi_balance.append(np.average(phi))
    for i in range(1, N-1):
        if phi[i] > max(phi[i-1], phi[i+1]):
            x_extrems.append(i)
            y_extrems.append(phi[i])
        
    # Создаем фигуру с тремя подграфиками
    fig, axs = plt.subplots(3, 1, figsize=(10, 12))
    
    # График энергии
    axs[0].plot(energies[:frame+1], 'r-')  # +1 потому что хотим включить текущий фрейм
    axs[0].set_xlim(0, frames-1)
    axs[0].set_ylim(min(energies)*0.9, max(energies)*1.1)
    axs[0].set_title("Энергия")

    # График баланса фазы
    axs[1].plot(phi_balance, 'r-') 
    axs[1].set_xlim(0, frames-1)
    axs[1].set_ylim(-0.5, 0.5)
    axs[1].set_title("Баланс фазы - avg(phi)")
    
    # График phi
    axs[2].plot(phi, 'r-')
    axs[2].plot(np.ones_like(phi), 'b--')
    axs[2].plot(-np.ones_like(phi), 'b--')
    axs[2].plot(x_extrems, y_extrems, 'bo')
    axs[2].set_xlim(0, N-1)
    axs[2].set_ylim(-1.2, 1.2)
    axs[2].set_title("Фaза")
    
    # Сохранение изображения
    # filename = f'gif_data/data{frame}.png'
    filename = f'gif_data/data_{frame}.png'
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    filenames.append(filename)

    print(f"{int(frame/frames*100)}%")
    x_extrems = []
    y_extrems = []

# Создание гифки
with imageio.get_writer(f'P-E_{eps}_{N}_{T_steps}.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

print(f"gif P-E_{eps}_{N}_{T_steps} had been created")

# Очистка временных файлов
# for filename in filenames:
#     os.remove(filename)