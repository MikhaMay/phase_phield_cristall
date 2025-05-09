import os
from dataclasses import dataclass
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import yaml


CUSTOM_RUN_DIR = None


@dataclass
class Parameters:
    time_steps: int
    output_interval: int
    # epsilon: float
    grid_size: float
    domain_length: float
    grid_spacing: float

    @classmethod
    def from_file(cls, path: Path) -> 'Parameters':
        with open(path, 'r') as file:
            params = yaml.safe_load(file)
        return cls(
            time_steps = params['timeSteps'],
            output_interval = params['outputInterval'],
            # epsilon = params['epsilon'],
            grid_size = params['gridSize'],
            domain_length = params['domainLength'],
            grid_spacing = params['gridSpacing'],
        )

latest_run_dir = CUSTOM_RUN_DIR or Path(
    Path('calculate_output_data/latest_run.txt').read_text().strip()
)
run_id = latest_run_dir.name
params = Parameters.from_file(latest_run_dir / 'params.yaml')


# Чтение энергий
energy_file = os.path.join(latest_run_dir, 'energies.bin')
with open(energy_file, 'rb') as file:
    N = np.fromfile(file, dtype=np.int32, count=1)[0]
    energies = np.fromfile(file, dtype=np.float64, count=N)

print(f'Loaded {len(energies)} energy values')

x_extrems = []
y_extrems = []
phi_balance = []

images_dir = Path('vizualize_output_data/images')
gifs_dir = Path('vizualize_output_data/gifs')

# Общее количество кадров
frames = (params.time_steps // params.output_interval) + 1
filenames = []
for frame in range(frames):
    data_file = os.path.join(latest_run_dir, f'{frame * params.output_interval}.bin')

    with open(data_file, 'rb') as file:
        N = np.fromfile(file, dtype=np.int32, count=1)[0]
        phi = np.fromfile(file, dtype=np.float64, count=N)
    phi_balance.append(np.average(phi))
    # Поиск экстремумов
    for i in range(1, N-1):
        if phi[i] > max(phi[i-1], phi[i+1]):
            x_extrems.append(i)
            y_extrems.append(phi[i])

    # Создаем фигуру с тремя подграфиками
    fig, axs = plt.subplots(3, 1, figsize=(10, 12))

    # График энергии
    axs[0].plot(energies[:frame+1], 'r-')
    axs[0].set_xlim(0, frames-1)
    if len(energies) > 0:
        axs[0].set_ylim(min(energies)*0.9, max(energies)*1.1)
    axs[0].set_title(f'Энергия')

    # График баланса фазы
    axs[1].plot(phi_balance, 'r-')
    axs[1].set_xlim(0, frames-1)
    axs[1].set_ylim(-0.5, 0.5)
    axs[1].set_title('Баланс фазы - avg(phi)')

    # График phi
    x_space = np.linspace(0, params.domain_length, N)
    axs[2].plot(x_space, phi, 'r-')
    axs[2].plot(x_space, np.ones_like(phi), 'b--')
    axs[2].plot(x_space, -np.ones_like(phi), 'b--')

    # Добавляем экстремумы, если они есть
    if x_extrems and y_extrems:
        x_extrems_space = [x_space[i] for i in x_extrems]
        axs[2].plot(x_extrems_space, y_extrems, 'bo')

    axs[2].set_xlim(0, params.domain_length)
    axs[2].set_ylim(-1.2, 1.2)
    axs[2].set_title('Фаза')
    axs[2].set_xlabel('x')

    # Сохранение изображения
    filename = images_dir / f'{frame}.png'
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    filenames.append(filename)

    print(f'Processed frame {frame+1}/{frames} ({int(frame/frames*100)}%)')
    x_extrems = []
    y_extrems = []

# Создание гифки
tag = os.path.basename(latest_run_dir)
gif_path = gifs_dir / f'{run_id}.gif'
print(f'Creating GIF: {gif_path}')

with imageio.get_writer(gif_path, mode='I', fps=10) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print(f'GIF "{gif_path}" created successfully')

# for filename in filenames:
#     os.remove(filename)