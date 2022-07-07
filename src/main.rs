use ggez::event::{self, EventHandler};
use ggez::graphics::{self, Color};
use ggez::{Context, ContextBuilder, GameResult};
use glam;

mod sim;
use sim::simulation::FluidCube;


fn translate_point(draw_param: &graphics::DrawParam, (x, y): (f32, f32)) -> glam::Vec2 {
    return glam::Vec2::from(draw_param.src.point()) + glam::vec2(x, y);
}

fn main() {
    // Make a Context.
    let mut cube = FluidCube::new(10, 1.0, 1.0, 1.0);
    cube.add_density(1, 1, 1, 1.0);
    cube.add_velocity(1, 1, 1, 1.0, 1.0, 1.0);
    cube.step();
    println!("Hello, world!");

    let (mut ctx, event_loop) = ContextBuilder::new("FluidSim", "Fluid Simulation")
        .build()
        .expect("ayooo, could not create context!");

    // Create an instance of your event handler.
    // Usually, you should provide it with the Context object to
    // use when setting your game up.
    let the_sim = SimulationState::new(&mut ctx);

    // Run!
    event::run(ctx, event_loop, the_sim);
}

struct SimulationState {
    // Your state here...
}

impl SimulationState {
    pub fn new(_ctx: &mut Context) -> SimulationState {
        // Load/create resources such as images here.
        SimulationState {
            // ...
        }
    }
}

impl EventHandler for SimulationState {
    fn update(&mut self, _ctx: &mut Context) -> GameResult<()> {
        // Update code here...
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