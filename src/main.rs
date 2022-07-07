// https://docs.rs/good-web-game/latest/good_web_game/

use good_web_game as ggez;
use std::path;
use std::env;

use ggez::event::{EventHandler};
use ggez::graphics::{self, Color};
use ggez::{Context, GameResult};
use glam;

mod sim;
use sim::simulation::FluidCube;


fn translate_point(draw_param: &graphics::DrawParam, (x, y): (f32, f32)) -> glam::Vec2 {
    return glam::Vec2::from(draw_param.src.point()) + glam::vec2(x, y);
}

fn main() -> GameResult {
    let mut cube = FluidCube::new(10, 1.0, 1.0, 1.0);
    cube.add_density(1, 1, 1, 1.0);
    cube.add_velocity(1, 1, 1, 1.0, 1.0, 1.0);
    cube.step();
    println!("Hello, world!");

    let resource_dir = if let Ok(manifest_dir) = env::var("CARGO_MANIFEST_DIR") {
        let mut path = path::PathBuf::from(manifest_dir);
        path.push("resources");
        path
    } else {
        path::PathBuf::from("./resources")
    };

    return ggez::start(
        ggez::conf::Conf::default()
        .high_dpi(true)
            // .cache(miniquad::conf::Cache::Tar(include_bytes!("resources.tar")))
            .physical_root_dir(Some(resource_dir)),
        |mut context| Box::new(SimulationState::new(&mut context)),
    );

}

struct SimulationState {
    
}

impl SimulationState {
    pub fn new(_ctx: &mut Context) -> SimulationState {
        SimulationState {
            
        }
    }
}

impl EventHandler for SimulationState {
    fn update(&mut self, _ctx: &mut Context) -> GameResult<()> {
        return Ok(());
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult<()> {
        graphics::clear(ctx, Color::BLACK);
        let (screen_width, screen_height) = graphics::drawable_size(ctx);
        let (screen_width_half, screen_height_half) = (screen_width * 0.5, screen_height * 0.5);
        let mut draw_param = graphics::DrawParam::default();

        // draw
        let text = graphics::Text::new("Hello, World!");
        let text_rect = text.dimensions(ctx);

        let text_pos = translate_point(
            &draw_param,
            (
                screen_width_half - (text_rect.w * 0.5),
                screen_height_half - (text_rect.h * 0.5)
            )
        );

        draw_param = draw_param.dest(text_pos);

        let _draw_result = graphics::draw(ctx, &text, draw_param);
        let _present_result = graphics::present(ctx);
        
        return Ok(());
    }
}
