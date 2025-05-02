import os
from dataclasses import dataclass
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib import cm


CUSTOM_RUN_DIR = None


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


latest_run_dir = CUSTOM_RUN_DIR or Path(
    Path('calculate_output_data/latest_run.txt').read_text().strip()
)
run_id = latest_run_dir.name
params = Parameters.from_file(latest_run_dir / 'params.yaml')


# Read energies
energy_file = os.path.join(latest_run_dir, 'energies.bin')
try:
    with open(energy_file, 'rb') as file:
        N = np.fromfile(file, dtype=np.int32, count=1)[0]
        # Second dimension is always 1 for energies
        M = np.fromfile(file, dtype=np.int32, count=1)[0]
        energies = np.fromfile(file, dtype=np.float64, count=N*M)
except FileNotFoundError as e:
    print(e)
    energies = []

print(f'Loaded {len(energies)} energy values')

phi_balance = []

images_dir = Path('vizualize_output_data/images')
images_dir.mkdir(parents=True, exist_ok=True)
gifs_dir = Path('vizualize_output_data/gifs')
gifs_dir.mkdir(parents=True, exist_ok=True)

# Total number of frames
frames = (params.time_steps // params.output_interval) + 1
filenames = []

# Create a regular grid for plotting
x = np.linspace(0, params.domain_length_x, params.grid_size_x)
y = np.linspace(0, params.domain_length_y, params.grid_size_y)
X, Y = np.meshgrid(x, y)

for frame in range(frames):
    data_file = os.path.join(latest_run_dir, f'{frame * params.output_interval}.bin')

    with open(data_file, 'rb') as file:
        nx = np.fromfile(file, dtype=np.int32, count=1)[0]
        ny = np.fromfile(file, dtype=np.int32, count=1)[0]
        
        # Read the 2D array data
        phi_raw = np.fromfile(file, dtype=np.float64, count=nx*ny)
        phi = phi_raw.reshape((nx, ny))

    # Calculate average phi for balance tracking
    phi_balance.append(np.mean(phi))

    # Create a figure with three subplots
    fig = plt.figure(figsize=(15, 10))
    
    # Top row: Energy and Phase Balance
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    # Bottom row: 2D phase plot
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    # Energy plot
    ax1.plot(energies[:frame+1], 'r-')
    ax1.set_xlim(0, frames-1)
    if len(energies) > 0:
        ax1.set_ylim(min(energies)*0.9, max(energies)*1.1)
    ax1.set_title('Energy')
    ax1.set_xlabel('Frame')
    ax1.set_ylabel('Energy')

    # Phase balance plot
    ax2.plot(phi_balance, 'r-')
    ax2.set_xlim(0, frames-1)
    ax2.set_ylim(-0.7, 0.1)
    ax2.set_title('Phase Balance - avg(phi)')
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Average φ')

    # 2D phase plot
    im = ax3.pcolormesh(X, Y, phi.T, cmap=cm.coolwarm, vmin=-1.0, vmax=1.0)
    ax3.set_title('Phase Field (φ)')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_aspect('equal')
    
    # Add a colorbar
    fig.colorbar(im, ax=ax3, label='φ value')

    # Add time information
    time = frame * params.output_interval * params.time_step
    fig.suptitle(f'Time: {time:.2f}, Frame: {frame}/{frames-1}', fontsize=16)

    # Tight layout with space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save the image
    filename = images_dir / f'{frame:04d}.png'
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    filenames.append(filename)

    print(f'Processed frame {frame+1}/{frames} ({int(frame/frames*100)}%)\t{filename}')

# Create GIF
tag = os.path.basename(latest_run_dir)
gif_path = gifs_dir / f'{run_id}.gif'
print(f'Creating GIF: {gif_path}')

with imageio.get_writer(gif_path, mode='I', fps=8) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print(f'GIF "{gif_path}" created successfully')

# Uncomment to delete image files after creating GIF
# for filename in filenames:
#     os.remove(filename)
