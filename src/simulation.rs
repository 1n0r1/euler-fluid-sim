use crate::cell::Cell;
use crate::cell::CellType;
use crate::space_domain::SpaceDomain;

use crate::presets;

const OMEGA: f32 = 1.7; // 0 <= OMEGA <= 2
const ITR_MAX: usize = 100;
const POISSON_EPSILON: f32 = 0.001;

pub struct Simulation {
    space_domain: SpaceDomain,

    delta_time: f32,        // seconds,
    acceleration: [f32; 2], // meters/seconds^2
    reynolds: f32,
    time: f32, // seconds
    initial_pressure_norm: Option<f32>,
    fluid_cell_count: Option<u32>,
}

impl Default for Simulation {
    fn default() -> Self {
        let preset = presets::cylinder_cross_flow();
        Simulation::from_preset(preset)
    }
}

impl Simulation {
    pub fn new() -> Self {
        let preset = presets::cylinder_cross_flow();
        Simulation::from_preset(preset)
    }

    pub fn from_preset(preset: presets::SimulationPreset) -> Self {
        Self {
            space_domain: preset.space_domain,
            delta_time: preset.delta_time,
            reynolds: preset.reynolds,
            acceleration: preset.acceleration,
            time: 0.0,
            initial_pressure_norm: None,
            fluid_cell_count: None,
        }
    }

    pub fn delta_space(&self) -> [f32; 2] {
        self.space_domain.delta_space()
    }

    pub fn space_size(&self) -> [usize; 2] {
        self.space_domain.space_size()
    }

    pub fn time(&self) -> f32 {
        self.time
    }

    pub fn pressure_range(&self) -> [f32; 2] {
        self.space_domain.pressure_range()
    }

    pub fn speed_range(&self) -> [f32; 2] {
        self.space_domain.speed_range()
    }

    pub fn psi_range(&self) -> [f32; 2] {
        self.space_domain.psi_range()
    }

    pub fn get_cell(&self, x: usize, y: usize) -> &Cell {
        self.space_domain.get_cell(x, y)
    }

    pub fn get_centered_velocity(&self, x: usize, y: usize) -> [f32; 2] {
        self.space_domain.get_centered_velocity(x, y)
    }

    pub fn iterate_one_timestep(&mut self) {
        // Change boundary cells and fluid cells next to boundary cells
        // velocity, pressure, f, g
        self.space_domain.update_boundary_velocities(); // O(n^2)
        self.space_domain.update_boundary_pressures_and_fg();

        // Change fluid cells f, g
        self.update_fg(); // O(n^2)

        // Change fluid cells rhs
        self.update_rhs(); // O(n^2)

        // Change fluid and boundary cells pressure
        self.solve_poisson_pressure_equation(); // O(m*n^2)

        // Change fluid cells velocity
        self.update_velocity(); // O(n^2)

        // Change psi of fluid cells and boundary cell on the left and bottom
        self.space_domain.update_psi(); // O(n^2)

        // For coloring
        self.space_domain.update_pressure_and_speed_range(); // O(n^2)

        self.time += self.delta_time
    }
}

impl Simulation {
    fn update_velocity(&mut self) {
        let space_size = self.space_domain.space_size();
        let delta_space = self.space_domain.delta_space();

        for x in 0..space_size[0] {
            for y in 0..space_size[1] {
                if let CellType::FluidCell = self.space_domain.get_cell(x, y).cell_type {
                    let right_cell_type: Option<CellType> = (x + 1 < space_size[0])
                        .then(|| self.space_domain.get_cell(x + 1, y).cell_type);
                    let top_cell_type: Option<CellType> = (y + 1 < space_size[1])
                        .then(|| self.space_domain.get_cell(x, y + 1).cell_type);

                    if let Some(right_cell_type) = right_cell_type {
                        if let CellType::BoundaryConditionCell(_) = right_cell_type {
                        } else {
                            self.space_domain.get_cell_mut(x, y).velocity[0] =
                                self.space_domain.get_cell(x, y).f
                                    - self.delta_time
                                        * (self.space_domain.get_cell(x + 1, y).pressure
                                            - self.space_domain.get_cell(x, y).pressure)
                                        / delta_space[0]
                        }
                    }

                    if let Some(top_cell_type) = top_cell_type {
                        if let CellType::BoundaryConditionCell(_) = top_cell_type {
                        } else {
                            self.space_domain.get_cell_mut(x, y).velocity[1] =
                                self.space_domain.get_cell(x, y).g
                                    - self.delta_time
                                        * (self.space_domain.get_cell(x, y + 1).pressure
                                            - self.space_domain.get_cell(x, y).pressure)
                                        / delta_space[1]
                        }
                    }
                }
            }
        }
    }

    fn get_initial_pressure_norm(&mut self) -> (f32, u32) {
        if let Some(x) = self.initial_pressure_norm {
            return (x, self.fluid_cell_count.unwrap());
        }
        let mut fluid_cell_count = 0;
        let space_size = self.space_domain.space_size();
        let mut initial_pressure_norm: f32 = (0..space_size[0])
            .map(|x| -> f32 {
                (0..space_size[1])
                    .map(|y| -> f32 {
                        match self.space_domain.get_cell(x, y).cell_type {
                            CellType::FluidCell => {
                                fluid_cell_count += 1;
                                (self.space_domain.get_cell(x, y)).pressure.powi(2)
                            }
                            _ => 0.0,
                        }
                    })
                    .sum()
            })
            .sum::<f32>()
            / (fluid_cell_count as f32);

        initial_pressure_norm = initial_pressure_norm.sqrt();
        self.initial_pressure_norm = Some(initial_pressure_norm);
        self.fluid_cell_count = Some(fluid_cell_count);
        (initial_pressure_norm, fluid_cell_count)
    }

    fn solve_poisson_pressure_equation(&mut self) {
        let space_size = self.space_domain.space_size();
        let delta_space = self.space_domain.delta_space();

        let (initial_pressure_norm, fluid_cell_count) = self.get_initial_pressure_norm();

        for _ in 0..ITR_MAX {
            let mut residual_norm: f32 = (0..space_size[0])
                .map(|x| -> f32 {
                    (0..space_size[1])
                        .map(|y| -> f32 {
                            if let CellType::FluidCell = self.space_domain.get_cell(x, y).cell_type
                            {
                                ((self.space_domain.get_cell(x + 1, y).pressure
                                    - 2.0 * self.space_domain.get_cell(x, y).pressure
                                    + self.space_domain.get_cell(x - 1, y).pressure)
                                    / delta_space[0].powi(2)
                                    + (self.space_domain.get_cell(x, y + 1).pressure
                                        - 2.0 * self.space_domain.get_cell(x, y).pressure
                                        + self.space_domain.get_cell(x, y - 1).pressure)
                                        / delta_space[1].powi(2)
                                    - self.space_domain.get_cell(x, y).rhs)
                                    .powi(2)
                            } else {
                                0.0
                            }
                        })
                        .sum()
                })
                .sum();

            residual_norm = (residual_norm / (fluid_cell_count as f32)).sqrt();

            if residual_norm < POISSON_EPSILON
                || residual_norm < initial_pressure_norm * POISSON_EPSILON
            {
                break;
            }

            self.update_pressures_for_boundary_cells();

            for x in 0..space_size[0] {
                for y in 0..space_size[1] {
                    if let CellType::FluidCell = self.space_domain.get_cell(x, y).cell_type {
                        self.space_domain.get_cell_mut(x, y).pressure = (1.0 - OMEGA)
                            * self.space_domain.get_cell(x, y).pressure
                            + OMEGA
                                * ((self.space_domain.get_cell(x + 1, y).pressure
                                    + (self.space_domain.get_cell(x - 1, y).pressure))
                                    / delta_space[0].powi(2)
                                    + (self.space_domain.get_cell(x, y + 1).pressure
                                        + (self.space_domain.get_cell(x, y - 1).pressure))
                                        / delta_space[1].powi(2)
                                    - self.space_domain.get_cell(x, y).rhs)
                                / (2.0 / delta_space[0].powi(2) + 2.0 / delta_space[1].powi(2));
                    }
                }
            }
        }
    }

    fn update_pressures_for_boundary_cells(&mut self) {
        let space_size = self.space_domain.space_size();

        for x in 0..space_size[0] {
            for y in 0..space_size[1] {
                let cell_type = &self.space_domain.get_cell(x, y).cell_type;

                if let CellType::BoundaryConditionCell(_) = cell_type {
                    let neighboring_cells = [
                        (x.wrapping_sub(1), y),
                        (x + 1, y),
                        (x, y.wrapping_sub(1)),
                        (x, y + 1),
                    ];

                    let mut neighboring_fluid_count = 0;
                    self.space_domain.get_cell_mut(x, y).pressure = 0.0;

                    for (dx, dy) in neighboring_cells.iter() {
                        if let Some(cell) = self.space_domain.try_get_cell(*dx, *dy) {
                            match cell.cell_type {
                                CellType::FluidCell => {
                                    self.space_domain.get_cell_mut(x, y).pressure += cell.pressure;
                                    neighboring_fluid_count += 1;
                                }
                                _ => {}
                            }
                        }
                    }

                    if neighboring_fluid_count != 0 {
                        self.space_domain.get_cell_mut(x, y).pressure =
                            self.space_domain.get_cell(x, y).pressure
                                / (neighboring_fluid_count as f32);
                    }
                }
            }
        }
    }

    fn update_rhs(&mut self) {
        let space_size = self.space_domain.space_size();
        let delta_space = self.space_domain.delta_space();

        for x in 0..space_size[0] {
            for y in 0..space_size[1] {
                if let CellType::FluidCell = self.space_domain.get_cell(x, y).cell_type {
                    self.space_domain.get_cell_mut(x, y).rhs =
                        ((self.space_domain.get_cell(x, y).f
                            - self.space_domain.get_cell(x - 1, y).f)
                            / delta_space[0]
                            + (self.space_domain.get_cell(x, y).g
                                - self.space_domain.get_cell(x, y - 1).g)
                                / delta_space[1])
                            / self.delta_time;
                }
            }
        }
    }

    fn update_fg(&mut self) {
        let space_size = self.space_domain.space_size();
        for x in 0..space_size[0] {
            for y in 0..space_size[1] {
                if let CellType::FluidCell = self.space_domain.get_cell(x, y).cell_type {
                    if let Some(CellType::FluidCell) = self
                        .space_domain
                        .try_get_cell(x + 1, y)
                        .map(|cell| cell.cell_type)
                    {
                        self.space_domain.get_cell_mut(x, y).f =
                            self.space_domain.get_cell(x, y).velocity[0]
                                + self.delta_time
                                    * ((self.space_domain.d2udx2(x, y)
                                        + self.space_domain.d2udy2(x, y))
                                        / self.reynolds
                                        - self.space_domain.du2dx(x, y)
                                        - self.space_domain.duvdy(x, y)
                                        + self.acceleration[0]);
                    }

                    if let Some(CellType::FluidCell) = self
                        .space_domain
                        .try_get_cell(x, y + 1)
                        .map(|cell| cell.cell_type)
                    {
                        self.space_domain.get_cell_mut(x, y).g =
                            self.space_domain.get_cell(x, y).velocity[1]
                                + self.delta_time
                                    * ((self.space_domain.d2vdx2(x, y)
                                        + self.space_domain.d2vdy2(x, y))
                                        / self.reynolds
                                        - self.space_domain.duvdx(x, y)
                                        - self.space_domain.dv2dy(x, y)
                                        + self.acceleration[1])
                    }
                }
            }
        }
    }
}
