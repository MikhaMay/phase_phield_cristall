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
            grid_size = params['gridSize'],
            domain_length = params['domainLength'],
            grid_spacing = params['gridSpacing'],
        )

latest_run_dir = CUSTOM_RUN_DIR or Path(
    Path('calculate_output_data/latest_run.txt').read_text().strip()
)
run_id = latest_run_dir.name
params = Parameters.from_file(latest_run_dir / 'params.yaml')


try:
    # Чтение энергий
    energy_file = os.path.join(latest_run_dir, 'energies.bin')
    with open(energy_file, 'rb') as file:
        N = np.fromfile(file, dtype=np.int32, count=1)[0]
        energies = np.fromfile(file, dtype=np.float64, count=N)

    print(f'Loaded {len(energies)} energy values')
except FileNotFoundError as e:
    energies = []
    print(str(e))

x_extrems = []
y_extrems = []
phi_balance = []
velocity_balance = []

images_dir = Path('vizualize_output_data/images')
gifs_dir = Path('vizualize_output_data/gifs')

# # Создаем директории, если они не существуют
# images_dir.mkdir(parents=True, exist_ok=True)
# gifs_dir.mkdir(parents=True, exist_ok=True)

dt = 1e-6

# Общее количество кадров
frames = (params.time_steps // params.output_interval) + 1
filenames = []
for frame in range(frames):
    step = frame * params.output_interval
    t = step * dt

    data_file = os.path.join(latest_run_dir, f'v_{step}.bin')

    # Чтение поля фазы
    with open(data_file, 'rb') as file:
        N = np.fromfile(file, dtype=np.int32, count=1)[0]
        phi = np.fromfile(file, dtype=np.float64, count=N)

    phi_balance.append(np.average(phi))

    # ОДИН график: фаза
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))

    x_space = np.linspace(0, params.domain_length, N)
    ax.plot(x_space, phi, 'r-', label='phi', linewidth=2)

    # БЕЗ экстремумов — ничего больше не рисуем поверх

    ax.set_xlim(0, params.domain_length)
    ax.set_ylim(4, 6)

    # Крупные и жирные надписи
    ax.set_title(rf'Скорость, $t = {t}$', fontsize=24, fontweight='bold')
    ax.set_xlabel('x', fontsize=20, fontweight='bold')
    ax.set_ylabel(r'$v$', fontsize=20, fontweight='bold')

    # Крупные подписи делений
    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()

    filename = images_dir / f'{frame}.png'
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    filenames.append(filename)

    print(f'Processed frame {frame+1}/{frames} ({int(frame/frames*100)}%)\t{filename}')


    # t_energy = np.arange(len(energies)) * params.output_interval * dt

    # fig, ax = plt.subplots(1, 1, figsize=(10, 4))

    # ax.plot(t_energy, energies, 'r-', linewidth=2)

    # ax.set_title(r'Энергия во времени', fontsize=24, fontweight='bold')
    # ax.set_xlabel(r'$t$', fontsize=20, fontweight='bold')
    # ax.set_ylabel(r'$E$', fontsize=20, fontweight='bold')

    # ax.tick_params(axis='both', which='major', labelsize=16)

    # plt.tight_layout()

    # energy_filename = images_dir / f'{run_id}_energy.png'
    # plt.savefig(energy_filename, dpi=150)
    # plt.close(fig)
