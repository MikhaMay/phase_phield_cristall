import argparse
import datetime
import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib import cm


@dataclass
class Parameters:
    time_step: float
    time_steps: int
    output_interval: int
    epsilon: float
    grid_size_x: int
    grid_size_y: int
    domain_length_x: float
    domain_length_y: float
    grid_spacing_x: float
    grid_spacing_y: float

    @classmethod
    def from_file(cls, path: Path) -> 'Parameters':
        with open(path, 'r') as file:
            params = yaml.safe_load(file)
        return cls(
            time_step=float(params['timeStep']),
            time_steps=params['timeSteps'],
            output_interval=params['outputInterval'],
            epsilon=params['epsilon'],
            grid_size_x=params['gridSizeX'],
            grid_size_y=params['gridSizeY'],
            domain_length_x=params['domainLengthX'],
            domain_length_y=params['domainLengthY'],
            grid_spacing_x=params['gridSpacingX'],
            grid_spacing_y=params['gridSpacingY'],
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
            M = np.fromfile(file, dtype=np.int32, count=1)[0]
            energies = np.fromfile(file, dtype=np.float64, count=N*M)
    except FileNotFoundError as e:
        print(e)
        energies = []
    
    print(f'Loaded {len(energies)} energy values')
    
    # Вычисляем общее количество фреймов
    frames = (params.time_steps // params.output_interval) + 1
    
    # Создаем массив для хранения баланса фазы
    phi_balance = []
    
    # Создаем сетку для отображения
    x = np.linspace(0, params.domain_length_x, params.grid_size_x)
    y = np.linspace(0, params.domain_length_y, params.grid_size_y)
    X, Y = np.meshgrid(x, y)
    
    # Обрабатываем все фреймы для расчета баланса фазы
    for frame in range(frames):
        data_file = simulation_dir / f'{frame * params.output_interval}.bin'
        
        # Если файл не существует, пропускаем
        if not data_file.exists():
            continue
        
        with open(data_file, 'rb') as file:
            nx = np.fromfile(file, dtype=np.int32, count=1)[0]
            ny = np.fromfile(file, dtype=np.int32, count=1)[0]
            
            # Читаем данные массива
            phi_raw = np.fromfile(file, dtype=np.float64, count=nx*ny)
            phi = phi_raw.reshape((nx, ny))
        
        # Рассчитываем средний баланс фазы
        phi_balance.append(np.mean(phi))
        
        # Если этот фрейм в списке выбранных, создаем визуализацию
        if frame in selected_frames:
            # Создаем и сохраняем визуализацию фазы
            plt.figure(figsize=(10, 8))
            im = plt.pcolormesh(X, Y, phi.T, cmap=cm.coolwarm, vmin=-1.0, vmax=1.0, shading='gouraud')
            plt.colorbar(im, label='φ value')
            plt.title(f'Phase Field (φ) - Frame {frame} - Time {frame * params.output_interval * params.time_step:.2f}')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
            plt.tight_layout()
            plt.savefig(output_dir / f'phase_frame_{frame:04d}.png', dpi=150)
            plt.close()
            
            print(f'Saved phase visualization for frame {frame}')
    
    # Создаем и сохраняем график энергии
    plt.figure(figsize=(12, 6))
    plt.plot(energies, 'r-')
    plt.title('Energy Plot')
    plt.xlabel('Frame')
    plt.ylabel('Energy')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / 'energy_plot.png', dpi=150)
    plt.close()
    
    print('Saved energy plot')
    
    # Создаем и сохраняем график баланса фазы
    plt.figure(figsize=(12, 6))
    plt.ylim(-0.7, 0.1)
    plt.plot(phi_balance, 'b-')
    plt.title('Phase Balance - avg(phi)')
    plt.xlabel('Frame')
    plt.ylabel('Average φ')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / 'phase_balance_plot.png', dpi=150)
    plt.close()
    
    print('Saved phase balance plot')
    
    print(f'All visualizations saved to {output_dir}')


simulation_dir = 'calculate_output_data/Lx20.0_Ly20.0_Nx64_Ny64_eps0.40_dt1.0e-08_Periodic_2d_random_-0.3avg'
process_simulation_data(simulation_dir, selected_frames=[0,10,50,100])
