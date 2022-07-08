// https://docs.rs/good-web-game/latest/good_web_game/

use good_web_game as ggez;

use ggez::event::EventHandler;
use ggez::graphics::{self, Color, DrawMode};
use ggez::{Context, GameResult};
use glam;

mod sim_2d;
use sim_2d::simulation;

use simulation::FluidPlane;
use simulation::N;
use simulation::SCALE;
use simulation::WINDOW_DIM;

fn main() -> GameResult {
    let resource_dir = if let Ok(manifest_dir) = std::env::var("CARGO_MANIFEST_DIR") {
        let mut path = std::path::PathBuf::from(manifest_dir);
        path.push("resources");
        path
    } else {
        std::path::PathBuf::from("./resources")
    };

    return ggez::start(
        ggez::conf::Conf::default()
            .window_width(WINDOW_DIM as i32)
            .window_height(WINDOW_DIM as i32)
            .window_resizable(false)
            .window_title("Fluid Simulation :D".to_owned())
            .high_dpi(true)
            // .cache(miniquad::conf::Cache::Tar(include_bytes!("resources.tar")))
            .physical_root_dir(Some(resource_dir)),
        |mut context| Box::new(SimulationState::new(&mut context)),
    );
}

struct SimulationState {
    size: i32,
    iter: i32,
    diffusion: f32,
    viscosity: f32,
    dt: f32,
    fluid: FluidPlane,
    scale: i32,
    real_size: i32,
    start_simulation: bool,
}

impl SimulationState {
    pub fn new(_ctx: &mut Context) -> SimulationState {
        let size = WINDOW_DIM;
        let diffusion = 0.0;
        let viscosity = 0.0;
        let dt = 0.1;
        let iter = 4;
        let fluid = FluidPlane::new(size as i32, SCALE as f32, diffusion, viscosity, dt, iter);
        SimulationState {
            size: size as i32,
            diffusion,
            viscosity,
            dt,
            fluid,
            iter,
            start_simulation: false,
            real_size: N as i32,
            scale: SCALE as i32,
        }
    }
}

impl EventHandler for SimulationState {
    fn mouse_motion_event(&mut self, _ctx: &mut Context, _x: f32, _y: f32, _dx: f32, _dy: f32) {
        self.start_simulation = true;
        let amount_x = _x - _dx;
        let amount_y = _y - _dy;
        self.fluid.add_density(_x as i32, _y as i32, 100.0);
        self.fluid.add_velocity(_x as i32 / self.scale, _y as i32 / self.scale, amount_x, amount_y);
        println!("\t {} {} {} {}", _x, _y, _dx, _dy)
    }

    fn update(&mut self, _ctx: &mut Context) -> GameResult<()> {
        return Ok(());
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult<()> {
        let canvas = graphics::Canvas::new(ctx, WINDOW_DIM as u16, WINDOW_DIM as u16)?;
        graphics::set_canvas(ctx, Some(&canvas));
        graphics::clear(ctx, Color::BLACK);

        if self.start_simulation {
            self.fluid.step();
            self.fluid.render(ctx, self.scale as f32);
            self.fluid.fade();
        }

        graphics::set_canvas(ctx, None);
        graphics::draw(ctx, &canvas, graphics::DrawParam::default())?;
        graphics::present(ctx)?;
        return Ok(());
    }
}
