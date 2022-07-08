
macro_rules! _IX {
    ($x:expr, $y:expr, $size:expr) => {
        (($x) + ($y) * ($size - 2)) as usize
    };
}

pub mod simulation {
    use good_web_game as ggez;
    use ggez::event::{EventHandler};
    use ggez::graphics::{self, Color, DrawMode};
    use ggez::{Context};
    use glam;

    use std::collections::HashMap;


    fn translate_point(draw_param: &graphics::DrawParam, (x, y): (f32, f32)) -> glam::Vec2 {
        return glam::Vec2::from(draw_param.src.point()) + glam::vec2(x, y);
    }

    fn clamp(v: f32, min: f32, max: f32) -> f32 {
        if v < min {
            return min;
        } else if v > max {
            return max;
        }

        return v;
    }

    // #[derive(Debug, Clone)]
    pub struct FluidPlane {
        iter: i32,
        size: i32,
        dt: f32,
        diff: f32,
        visc: f32,
        s: Vec<f32>,
        density: Vec<f32>,
        vx: Vec<f32>,
        vy: Vec<f32>,
        vx0: Vec<f32>,
        vy0: Vec<f32>,

        rect_cache: HashMap<(i32, i32), graphics::Rect>,
    }

    impl FluidPlane {
        pub fn new(size: i32, scale: f32, diffusion: f32, viscosity: f32, dt: f32, iter: i32) -> Self {
            let arr_size: usize = size as usize * size as usize;
            let mut rect_cache = HashMap::new();

            for i in 0..size {
                for j in 0..size {
                    let x = i as f32 * scale;
                    let y = j as f32 * scale;
                    let rect = graphics::Rect::new(x, y, scale, scale);

                    rect_cache.insert((i, j), rect);

                }
            }

            return FluidPlane {
                size: size,
                dt: dt,
                diff: diffusion,
                visc: viscosity,
                s: vec![0_f32; arr_size],
                density: vec![0_f32; arr_size],
                vx: vec![0_f32; arr_size],
                vy: vec![0_f32; arr_size],
                vx0: vec![0_f32; arr_size],
                vy0: vec![0_f32; arr_size],
                iter: iter,
                rect_cache,
            };
        }

        pub fn step(&mut self) {
            let size_n = self.size;
            let visc = self.visc;
            let diff = self.diff;
            let dt = self.dt;
            let iter = self.iter;
            let mut vx = &mut self.vx;
            let mut vy = &mut self.vy;
            let mut vx0 = &mut self.vx0;
            let mut vy0 = &mut self.vy0;
            let mut s = &mut self.s;
            let mut density = &mut self.density;
            utils::diffuse(1, &mut vx0, &vx, visc, dt, iter, size_n);
            utils::diffuse(2, &mut vy0, &vy, visc, dt, iter, size_n);
            utils::project(&mut vx0, &mut vy0, &mut vx, &mut vy, iter, size_n);
            utils::advect(1, &mut vx, &vx0, &vx0, &vy0, dt, size_n);
            utils::advect(2, &mut vy, &vy0, &vx0, &vy0, dt, size_n);
            utils::project(&mut vx, &mut vy, &mut vx0, &mut vy0, iter, size_n);
            utils::diffuse(0, &mut s, &density, diff, dt, iter, size_n);
            utils::advect(0, &mut density, &s, &vx, &vy, dt, size_n);
        }

        pub fn add_density(&mut self, x: i32, y: i32, amount: f32) {
            let size_n = self.size;
            self.density[_IX!(x, y, size_n)] += amount;
        }

        pub fn add_velocity(&mut self, x: i32, y: i32, amount_x: f32, amount_y: f32) {
            let size_n = self.size;
            let index = _IX!(x, y, size_n);
            self.vx[index] += amount_x;
            self.vy[index] += amount_y;
        }

        pub fn render(&self, ctx: &mut Context, scale_float: f32) {
            let mut draw_param = graphics::DrawParam::default();
            let draw_mode = DrawMode::Fill(graphics::FillOptions::DEFAULT);
            for i in 0..self.size {
                for j in 0..self.size {
                    let x = i as f32 * scale_float;
                    let y = j as f32 * scale_float;
                    let d = self.density[_IX!(i, j, self.size)];

                    draw_param = draw_param.dest(translate_point(&draw_param, (x, y)));
                    let rect = *self.rect_cache.get(&(i, j)).unwrap();
                    let mesh = graphics::MeshBuilder::new()
                        .rectangle(draw_mode, rect, Color::new(1.0, 1.0, 1.0, d))
                        .unwrap()
                        .build(ctx)
                        .unwrap();

                    graphics::draw(
                        ctx,
                        &mesh,
                        draw_param,
                    ).unwrap();
                }
            }
        }

        pub fn fade(&mut self) {
            for i in 0..self.density.len() {
                let d = self.density[i];
                self.density[i] = clamp(d - 0.1, 0.0, 1.0);
            }
        }
    }
    pub mod utils {
        fn set_bnd(b: i32, x: &mut Vec<f32>, size_n: i32) {
            for i in 1..(size_n - 1) {
                x[_IX!(i, 0, size_n)] = if b == 2 {
                    -x[_IX!(i, 1, size_n)]
                } else {
                    x[_IX!(i, 1, size_n)]
                };
                x[_IX!(i, size_n - 1, size_n)] = if b == 2 {
                    -x[_IX!(i, size_n - 2, size_n)]
                } else {
                    x[_IX!(i, size_n - 2, size_n)]
                };
            }

            for j in 1..(size_n - 1) {
                x[_IX!(0, j, size_n)] = if b == 1 {
                    -x[_IX!(1, j, size_n)]
                } else {
                    x[_IX!(1, j, size_n)]
                };
                x[_IX!(size_n - 1, j, size_n)] = if b == 1 {
                    -x[_IX!(size_n - 2, j, size_n)]
                } else {
                    x[_IX!(size_n - 2, j, size_n)]
                };
            }

            x[_IX!(0, 0, size_n)] = 0.5_f32 * (x[_IX!(1, 0, size_n)] + x[_IX!(0, 1, size_n)]);
            x[_IX!(0, size_n - 1, size_n)] =
                0.5_f32 * (x[_IX!(1, size_n - 1, size_n)] + x[_IX!(0, size_n - 2, size_n)]);
            x[_IX!(size_n - 1, 0, size_n)] =
                0.5_f32 * (x[_IX!(size_n - 2, 0, size_n)] + x[_IX!(size_n - 1, 1, size_n)]);
            x[_IX!(size_n - 1, size_n - 1, size_n)] = 0.5_f32
                * (x[_IX!(size_n - 2, size_n - 1, size_n)]
                    + x[_IX!(size_n - 1, size_n - 2, size_n)]);
        }

        fn lin_solve(
            b: i32,
            x: &mut Vec<f32>,
            x0: &Vec<f32>,
            a: f32,
            c: f32,
            iter: i32,
            size_n: i32,
        ) {
            let c_recip: f32 = 1.0_f32 / c;
            for _ in 0..iter {
                for j in 1..(size_n - 1) {
                    for i in 1..(size_n - 1) {
                        x[_IX!(i, j, size_n)] = (x0[_IX!(i, j, size_n)]
                            + a * (x[_IX!(i + 1, j, size_n)]
                                + x[_IX!(i - 1, j, size_n)]
                                + x[_IX!(i, j + 1, size_n)]
                                + x[_IX!(i, j - 1, size_n)]))
                            * c_recip;
                    }
                }
                set_bnd(b, x, size_n);
            }
        }

        pub fn diffuse(
            b: i32,
            x: &mut Vec<f32>,
            x0: &Vec<f32>,
            diff: f32,
            dt: f32,
            iter: i32,
            size_n: i32,
        ) {
            let size_n_f32: f32 = size_n as f32;
            let a: f32 = dt * diff * (size_n_f32 - 2.0) * (size_n_f32 - 2.0);
            lin_solve(b, x, x0, a, 1.0 + 6.0 * a, iter, size_n);
        }

        pub fn advect(
            b: i32,
            d: &mut Vec<f32>,
            d0: &Vec<f32>,
            veloc_x: &Vec<f32>,
            veloc_y: &Vec<f32>,
            dt: f32,
            size_n: i32,
        ) {
            let n_float: f32 = size_n as f32;
            let mut i0: f32;
            let mut i1: f32;
            let mut j0: f32;
            let mut j1: f32;

            let dtx: f32 = dt * (n_float - 2.0);
            let dty: f32 = dt * (n_float - 2.0);

            let mut s0: f32;
            let mut s1: f32;
            let mut t0: f32;
            let mut t1: f32;

            let mut tmp1: f32;
            let mut tmp2: f32;

            let mut x: f32;
            let mut y: f32;

            let mut ifloat: f32 = 1.0;
            let mut jfloat: f32 = 1.0;

            let mut i: i32 = 1;
            let mut j: i32 = 1;

            while j < (size_n - 1) {
                while i < (size_n - 1) {
                    tmp1 = dtx * veloc_x[_IX!(i, j, size_n)];
                    tmp2 = dty * veloc_y[_IX!(i, j, size_n)];
                    x = ifloat - tmp1;
                    y = jfloat - tmp2;
                    if x < 0.5_f32 {
                        x = 0.5_f32;
                    }
                    if x > n_float + 0.5_f32 {
                        x = n_float + 0.5_f32;
                    }
                    i0 = x.floor();
                    i1 = i0 + 1.0_f32;
                    if y < 0.5_f32 {
                        y = 0.5_f32;
                    }
                    if y > n_float + 0.5_f32 {
                        y = n_float + 0.5_f32;
                    }
                    j0 = y.floor();
                    j1 = j0 + 1.0_f32;

                    s1 = x - i0;
                    s0 = 1.0_f32 - s1;
                    t1 = y - j0;
                    t0 = 1.0_f32 - t1;

                    let i0i: i32 = i0 as i32;
                    let i1i: i32 = i1 as i32;
                    let j0i: i32 = j0 as i32;
                    let j1i: i32 = j1 as i32;

                    d[_IX!(i, j, size_n)] = 
                        s0 * (t0 * d0[_IX!(i0i, j0i, size_n)] + (t1 * (d0[_IX!(i0i, j1i, size_n)])))+ 
                        s1 * (t0 * d0[_IX!(i1i, j0i, size_n)] + (t1 * (d0[_IX!(i1i, j1i, size_n)])));
                    i += 1;
                    ifloat += 1.0;
                }
                j += 1;
                jfloat += 1.0;
            }
            set_bnd(b, d, size_n);
        }

        pub fn project(
            veloc_x: &mut Vec<f32>,
            veloc_y: &mut Vec<f32>,
            p: &mut Vec<f32>,
            div: &mut Vec<f32>,
            iter: i32,
            size_n: i32,
        ) {
            let float_n: f32 = size_n as f32;
            for j in 1..(size_n - 1) {
                for i in 1..(size_n - 1) {
                    div[_IX!(i, j, size_n)] = -0.5
                        * (veloc_x[_IX!(i + 1, j, size_n)] - veloc_x[_IX!(i - 1, j, size_n)]
                            + veloc_y[_IX!(i, j + 1, size_n)]
                            - veloc_y[_IX!(i, j - 1, size_n)])
                        / float_n;
                    p[_IX!(i, j, size_n)] = 0.0;
                }
            }
            set_bnd(0, div, size_n);
            set_bnd(0, p, size_n);
            lin_solve(0, p, div, 1.0, 6.0, iter, size_n);
            for j in 1..(size_n - 1) {
                for i in 1..(size_n - 1) {
                    veloc_x[_IX!(i, j, size_n)] -=
                        0.5 * (p[_IX!(i + 1, j, size_n)] - p[_IX!(i - 1, j, size_n)]) * float_n;
                    veloc_y[_IX!(i, j, size_n)] -=
                        0.5 * (p[_IX!(i, j + 1, size_n)] - p[_IX!(i, j - 1, size_n)]) * float_n;
                }
            }
            set_bnd(1, veloc_x, size_n);
            set_bnd(2, veloc_y, size_n);
        }
    }
}
