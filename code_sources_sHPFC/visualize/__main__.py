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
            time_steps=params['timeSteps'],
            output_interval=params['outputInterval'],
            grid_size=params['gridSize'],
            domain_length=params['domainLength'],
            grid_spacing=params['gridSpacing'],
        )


def read_binary_field(path: Path) -> np.ndarray:
    with open(path, 'rb') as file:
        n = np.fromfile(file, dtype=np.int32, count=1)[0]
        data = np.fromfile(file, dtype=np.float64, count=n)
    return data


latest_run_dir = CUSTOM_RUN_DIR or Path(
    Path('calculate_output_data/latest_run.txt').read_text().strip()
)
run_id = latest_run_dir.name
params = Parameters.from_file(latest_run_dir / 'params.yaml')


try:
    energy_file = latest_run_dir / 'energies.bin'
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

# Если нужно — раскомментируй
# images_dir.mkdir(parents=True, exist_ok=True)
# gifs_dir.mkdir(parents=True, exist_ok=True)

frames = (params.time_steps // params.output_interval) + 1
filenames = []

for frame in range(frames):
    step = frame * params.output_interval

    data_file = latest_run_dir / f'{step}.bin'
    xi_file = latest_run_dir / f'xi_{step}.bin'
    coarse_xi_file = latest_run_dir / f'coarse_xi_{step}.bin'
    v_file = latest_run_dir / f'v_{step}.bin'
    vlap_file = latest_run_dir / f'vlap_{step}.bin'
    f_ext_file = latest_run_dir / f'f_ext_{step}.bin'

    phi = read_binary_field(data_file)
    xi = read_binary_field(xi_file)
    coarse_xi = read_binary_field(coarse_xi_file)
    v = read_binary_field(v_file)
    vlap = read_binary_field(vlap_file)
    f_ext = read_binary_field(f_ext_file)

    N = len(phi)

    phi_balance.append(np.average(phi))
    velocity_balance.append(np.average(v))

    for i in range(1, N - 1):
        if phi[i] > max(phi[i - 1], phi[i + 1]):
            x_extrems.append(i)
            y_extrems.append(phi[i])

    # Теперь 8 графиков
    fig, axs = plt.subplots(8, 1, figsize=(10, 18))

    # Энергия
    axs[0].plot(energies[:frame + 1], 'r-')
    axs[0].set_xlim(0, frames - 1)
    if len(energies) > 0:
        axs[0].set_ylim(min(energies) * 0.9, max(energies) * 1.1)
    axs[0].set_title('Энергия')
    axs[0].set_xlabel('t')

    # Баланс фазы
    axs[1].plot(phi_balance, 'r-')
    axs[1].set_xlim(0, frames - 1)
    axs[1].set_ylim(-0.5, 0.5)
    axs[1].set_title('Баланс фазы - avg(phi)')
    axs[1].set_xlabel('t')

    x_space = np.linspace(0, params.domain_length, N)

    # phi
    axs[2].plot(x_space, phi, 'r-')
    axs[2].plot(x_space, np.ones_like(phi), 'b--')
    axs[2].plot(x_space, -np.ones_like(phi), 'b--')

    if x_extrems and y_extrems:
        x_extrems_space = [x_space[i] for i in x_extrems]
        axs[2].plot(x_extrems_space, y_extrems, 'bo')

    axs[2].set_xlim(0, params.domain_length)
    axs[2].set_ylim(-1.2, 1.2)
    axs[2].set_title('Фаза')
    axs[2].set_xlabel('x')

    # f_ext — прямо под phi
    axs[3].plot(x_space, f_ext, 'm-')
    axs[3].set_xlim(0, params.domain_length)

    f_ext_abs_max = np.max(np.abs(f_ext)) if len(f_ext) > 0 else 0.0
    if f_ext_abs_max > 0:
        axs[3].set_ylim(-1.1 * f_ext_abs_max, 1.1 * f_ext_abs_max)
    else:
        axs[3].set_ylim(-1.0, 1.0)

    axs[3].set_title('Внешняя сила $f_{ext}$')
    axs[3].set_xlabel('x')

    # v — под графиком силы
    axs[4].plot(x_space, v, 'r-')
    axs[4].set_xlim(0, params.domain_length)
    axs[4].set_ylim(-5, 5)
    axs[4].set_title('Поле скорости')
    axs[4].set_xlabel('x')

    # xi
    axs[5].plot(x_space, xi, 'r-')
    axs[5].set_xlim(0, params.domain_length)
    axs[5].set_ylim(-5, 5)
    axs[5].set_title('Поле xi')
    axs[5].set_xlabel('x')

    # coarse_xi
    axs[6].plot(x_space, coarse_xi, 'r-')
    axs[6].set_xlim(0, params.domain_length)
    axs[6].set_ylim(-0.2, 0.2)
    axs[6].set_title('Поле coarse_xi')
    axs[6].set_xlabel('x')

    # vlap
    axs[7].plot(x_space, vlap, 'r-')
    axs[7].set_xlim(0, params.domain_length)
    axs[7].set_ylim(-5, 5)
    axs[7].set_title('Поле лапласиана скорости')
    axs[7].set_xlabel('x')

    filename = images_dir / f'{frame}.png'
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)
    filenames.append(filename)

    print(f'Processed frame {frame + 1}/{frames} ({int(frame / frames * 100)}%)\t{filename}')
    x_extrems = []
    y_extrems = []

gif_path = gifs_dir / f'{run_id}.gif'
print(f'Creating GIF: {gif_path}')

with imageio.get_writer(gif_path, mode='I', fps=10) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print(f'GIF "{gif_path}" created successfully')
