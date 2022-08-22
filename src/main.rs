// https://github.com/bevyengine/bevy/tree/latest/examples/2d
// https://bevy-cheatbook.github.io/platforms/wasm.html

mod utils;
mod fluids;
mod debug;

use fluids::simulation::sim_2d::{FluidPlane, index, clamp_index};
use debug::DebugPlugin;

use bevy::prelude::*;
use bevy::input::mouse::MouseMotion;
use bevy::sprite::MaterialMesh2dBundle;
use bevy_inspector_egui::Inspectable;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugin(DebugPlugin)
        .add_startup_system_to_stage(StartupStage::PostStartup, spawn_fluid_system)
        .add_system(user_mouse_fluid_system)
        .add_system(simulate_fluid_system)
        .add_system(display_fluid_state_system)
        .add_system(fade_fluid_system)
        .run();
}


fn spawn_fluid_system(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    let sim_state = SimulationState::default();

    commands.spawn_bundle(Camera2dBundle::default());
    let tile_size = 10.0;
    let sz = sim_state.size;
    for i in 0..sz {
        for j in 0..sz {
            let x = i as f32;
            let y = j as f32;

            let mesh = MaterialMesh2dBundle {
                mesh: meshes.add(Mesh::from(shape::Quad::default())).into(),
                transform: Transform::from_xyz(x * tile_size, -y * tile_size, 100.).with_scale(Vec3::splat(tile_size)),
                material: materials.add(ColorMaterial::from(Color::rgba(1., 1., 1., 1.))),
                ..default()
            };
            commands.spawn_bundle(mesh);
        }
    }

    commands.spawn().insert(sim_state);
}

fn user_mouse_fluid_system(
    // mut commands: Commands,
    windows: Res<Windows>,
    mut motion_evr: EventReader<MouseMotion>,
    mut cursor_evr: EventReader<CursorMoved>,
    mut query: Query<&mut SimulationState>,
) {
    let window = windows.get_primary().expect("No primary display");

    if let Some(_position) = window.cursor_position() {
        let movement_delta_opt = motion_evr.iter().last();
        let target_pos_opt = cursor_evr.iter().last();

        if let Some(cursor_dxy) = movement_delta_opt {
            if let Some(cursor_pos) = target_pos_opt {
                let (dx, dy) = (cursor_dxy.delta.x, cursor_dxy.delta.y);
                let (px, py) = (cursor_pos.position.x, cursor_pos.position.y);

                for mut sim_state in query.iter_mut() {
                    sim_state.start_simulation = true;
                    let scale = sim_state.scale;
                    sim_state.fluid.add_density(px as i32, py as i32, 100.0);
                    sim_state.fluid.add_velocity(
                        px as i32 / scale,
                        py as i32 / scale,
                        dx,
                        dy,
                    );
                }
            }
        }
    }
}

fn simulate_fluid_system(
    mut query: Query<&mut SimulationState>,
) {
    for mut sim_state in query.iter_mut() {
        sim_state.fluid.step();
    }
}

fn fade_fluid_system(
    mut query: Query<&mut SimulationState>,
) {
    for mut sim_state in query.iter_mut() {
        sim_state.fluid.fade();
    }
}

fn display_fluid_state_system(
    // time: Res<Time>,
    mut query: Query<&mut SimulationState>
) {
    for sim_state in query.iter_mut() {
        if sim_state.start_simulation {
            let sz = sim_state.size;
            for i in 0..sz {
                for j in 0..sz {
                    let _x = i as f32;
                    let _y = j as f32;
                    let _d = sim_state.fluid.density[index!(i, j, sz)];

                    // let p = graphics::DrawParam::new()
                    //     .color(graphics::Color::new(1.0, 1.0, 1.0, d * scale_float))
                    //     .dest(glam::Vec2::new(x, y));
                    // mesh_batch.add(p);
                }
            }
        }
    }
}

#[derive(Debug, Component, Inspectable)]
struct SimulationState {
    size: i32,
    iter: i32,
    diffusion: f32,
    viscosity: f32,
    dt: f32,
    scale: i32,
    real_size: i32,
    start_simulation: bool,
    pub fluid: FluidPlane,
}

impl Default for SimulationState {
    fn default() -> Self {
        let size = 42;
        let diffusion = 0.0;
        let viscosity = 0.0;
        let dt = 0.1;
        let iter = 4;
        let scale = 10.0;
        let fluid = FluidPlane::new(
            size as i32,
            scale as f32,
            diffusion,
            viscosity,
            dt,
            iter,
        );

        Self {
            size: size as i32,
            diffusion,
            viscosity,
            dt,
            fluid,
            iter,
            start_simulation: true,
            real_size: size as i32,
            scale: scale as i32,
        }
    }
}
