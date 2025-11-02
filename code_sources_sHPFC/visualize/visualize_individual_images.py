import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import yaml


LABELSIZE = 18
FONTSIZE = 24
TAG = 'v_'

@dataclass
class Parameters:
    time_step: float
    time_steps: int
    output_interval: int
    # epsilon: float
    grid_size: int
    domain_length: float
    grid_spacing: float

    @classmethod
    def from_file(cls, path: Path) -> 'Parameters':
        with open(path, 'r') as file:
            params = yaml.safe_load(file)
        return cls(
            time_step=float(params['timeStep']),
            time_steps=params['timeSteps'],
            output_interval=params['outputInterval'],
            # epsilon=params['epsilon'],
            grid_size=params['gridSize'],
            domain_length=params['domainLength'],
            grid_spacing=params['gridSpacing'],
        )

def process_simulation_data(simulation_dir_path, selected_frames):
    # Преобразуем путь в объект Path, если он передан как строка
    simulation_dir = Path(simulation_dir_path)
    run_id = simulation_dir.name

    # Загружаем параметры симуляции
    params = Parameters.from_file(simulation_dir / 'params.yaml')

    # Создаем уникальное имя для директории с результатами
    output_dir_name = f"images_{run_id}"
    output_dir = Path(f"separate_images/{output_dir_name}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Записываем информацию о запуске
    with open(output_dir / "info.txt", "w") as f:
        f.write(f"Source directory: {simulation_dir}\n")
        f.write(f"Run ID: {run_id}\n")
        f.write(f"Parameters: \n")
        for key, value in params.__dict__.items():
            f.write(f"  {key}: {value}\n")

    # Читаем данные энергии
    energy_file = simulation_dir / 'energies.bin'
    try:
        with open(energy_file, 'rb') as file:
            N = np.fromfile(file, dtype=np.int32, count=1)[0]
            energies = np.fromfile(file, dtype=np.float64, count=N)
    except FileNotFoundError as e:
        print(e)
        energies = []

    print(f'Loaded {len(energies)} energy values')

    # Вычисляем общее количество фреймов
    frames = (params.time_steps // params.output_interval) + 1
    print(frames)

    # Создаем массив для хранения баланса фазы
    phi_balance = []

    # Создаем сетку для отображения
    x = np.linspace(0, params.domain_length, params.grid_size)
    times = []

    # Обрабатываем все фреймы для расчета баланса фазы
    for frame in range(frames):
        times.append(frame * params.output_interval * params.time_step)
        data_file = simulation_dir / f'{TAG}{frame * params.output_interval}.bin'

        with open(data_file, 'rb') as file:
            nx = np.fromfile(file, dtype=np.int32, count=1)[0]

            # Читаем данные массива
            phi_raw = np.fromfile(file, dtype=np.float64, count=nx)
            phi = phi_raw

        # Рассчитываем средний баланс фазы
        phi_balance.append(np.mean(phi))

        # Если этот фрейм не в списке выбранных, пропускаем визуализацию
        if frame not in selected_frames:
            continue
        # Создаем и сохраняем визуализацию фазы
        plt.figure(figsize=(10, 8))
        plt.plot(x, phi, 'r-')
        plt.xlabel('x', fontsize=FONTSIZE)
        plt.ylabel('v value', fontsize=FONTSIZE)
        plt.ylim(-4.0, 4.0)
        plt.xlim(0.0, 80.0)
        plt.xticks(fontsize=LABELSIZE)
        plt.yticks(fontsize=LABELSIZE)
        plt.tight_layout()
        filename = output_dir / f'{TAG}phase_frame_{frame:04d}.png'
        plt.savefig(filename, dpi=300)
        plt.close()

        print(f'Saved phase visualization for frame {frame}')

        with imageio.get_writer(f'{str(filename)}.gif', mode='I', fps=1) as writer:
            image = imageio.imread(filename)
            writer.append_data(image)

        subprocess.run(['convert', f'{str(filename)}.gif', str(output_dir/f'frame_{frame}.png')])
        subprocess.run(['rm', f'{str(filename)}.gif'])
        print(f'converted {frame}')

    # Создаем и сохраняем график энергии
    plt.figure(figsize=(8, 6))
    plt.plot(times, energies, 'r-')
    plt.xlim(times[0], times[-1])
    plt.ylim(0, 20)
    plt.xticks(fontsize=LABELSIZE)
    plt.yticks(fontsize=LABELSIZE)
    plt.xlabel('Время', fontsize=FONTSIZE)
    plt.ylabel('Свободная энергия', fontsize=FONTSIZE)
    plt.tight_layout()
    plt.savefig(output_dir / 'energy_plot.png', dpi=150)
    plt.close()

    print('Saved energy plot')

    # Создаем и сохраняем график баланса фазы
    plt.figure(figsize=(8, 6))
    plt.ylim(-0.7, 0.1)
    plt.plot(times, phi_balance, 'r-')
    plt.xlim(times[0], times[-1])
    plt.xticks(fontsize=LABELSIZE)
    plt.yticks(fontsize=LABELSIZE)
    plt.xlabel('Время', fontsize=FONTSIZE)
    plt.ylabel('Общая фазовая концентрация', fontsize=FONTSIZE)
    plt.tight_layout()
    plt.savefig(output_dir / 'phase_balance_plot.png', dpi=150)
    plt.close()

    print('Saved phase balance plot')

    print(f'All visualizations saved to {output_dir}')

simulation_dir = 'calculate_output_data/L80.0_N200_r-0.40_dt1.0e-04_Periodic_sHPFC'
process_simulation_data(simulation_dir, selected_frames=[10, 30, 100])
