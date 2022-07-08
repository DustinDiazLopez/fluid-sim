macro_rules! _IX {
    ($x:expr, $y:expr, $size:expr) => {
        (($x) + ($y) * ($size)) as usize
    };
}

pub mod simulation {
    use ggez::event::EventHandler;
    use ggez::graphics::{self, Color, DrawMode};
    use ggez::Context;
    use glam;
    use good_web_game as ggez;
    use std::collections::HashMap;
    pub const N: usize = 32;
    pub const SCALE: usize = 10;
    pub const WINDOW_DIM: usize = N * SCALE;
    pub const ARRAY_SIZE: usize = WINDOW_DIM * WINDOW_DIM;

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

    fn clampi(v: i32, min: i32, max: i32) -> i32 {
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
        s: [f32; ARRAY_SIZE],
        density: [f32; ARRAY_SIZE],
        vx: [f32; ARRAY_SIZE],
        vy: [f32; ARRAY_SIZE],
        vx0: [f32; ARRAY_SIZE],
        vy0: [f32; ARRAY_SIZE],
    }

    impl FluidPlane {
        pub fn new(
            size: i32,
            scale: f32,
            diffusion: f32,
            viscosity: f32,
            dt: f32,
            iter: i32,
        ) -> Self {
            let arr_size: usize = size as usize * size as usize;

            return FluidPlane {
                size: size,
                dt: dt,
                diff: diffusion,
                visc: viscosity,
                s: [0_f32; ARRAY_SIZE],
                density: [0_f32; ARRAY_SIZE],
                vx: [0_f32; ARRAY_SIZE],
                vy: [0_f32; ARRAY_SIZE],
                vx0: [0_f32; ARRAY_SIZE],
                vy0: [0_f32; ARRAY_SIZE],
                iter: iter,
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
            diffuse(1, &mut vx0, &mut vx, visc, dt, iter, size_n);
            diffuse(2, &mut vy0, &mut vy, visc, dt, iter, size_n);
            project(&mut vx0, &mut vy0, &mut vx, &mut vy, iter, size_n);
            advect(1, &mut vx, &vx0, &vx0, &mut vy0, dt, size_n);
            advect(2, &mut vy, &vy0, &vx0, &vy0, dt, size_n);
            project(&mut vx, &mut vy, &mut vx0, &mut vy0, iter, size_n);
            diffuse(0, &mut s, &mut density, diff, dt, iter, size_n);
            advect(0, &mut density, &mut s, &mut vx, &mut vy, dt, size_n);
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

        pub fn render(&mut self, ctx: &mut Context, scale_float: f32) {
            let filled_rect = graphics::MeshBuilder::new()
                .rectangle(
                    graphics::DrawMode::fill(),
                    graphics::Rect::new(0.0, 0.0, scale_float, scale_float),
                    graphics::Color::WHITE,
                )
                .unwrap()
                .build(ctx)
                .unwrap();
            let mut mesh_batch = graphics::MeshBatch::new(filled_rect).unwrap();
            for i in 0..self.size {
                for j in 0..self.size {
                    let x = i as f32;
                    let y = j as f32;
                    let d = self.density[_IX!(i, j, self.size)];

                    // let rect = graphics::Rect::new(x, y, scale_float, scale_float);
                    let p = graphics::DrawParam::new()
                        .color(graphics::Color::new(1.0, 1.0, 1.0, d * scale_float))
                        .dest(glam::Vec2::new(x, y));
                    mesh_batch.add(p);
                }
            }

            mesh_batch
                .draw(ctx, graphics::DrawParam::default())
                .unwrap();
        }

        pub fn fade(&mut self) {
            for i in 0..self.density.len() {
                let d = self.density[i];
                self.density[i] = clamp(d - 0.1, 0.0, 1.0);
            }
        }
    }
    fn set_bnd(b: i32, x: &mut [f32], size_n: i32) {
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
            * (x[_IX!(size_n - 2, size_n - 1, size_n)] + x[_IX!(size_n - 1, size_n - 2, size_n)]);
    }

    fn lin_solve(b: i32, x: &mut [f32], x0: &[f32], a: f32, c: f32, iter: i32, size_n: i32) {
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
        x: &mut [f32; ARRAY_SIZE],
        x0: &[f32; ARRAY_SIZE],
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
        d: &mut [f32; ARRAY_SIZE],
        d0: &[f32; ARRAY_SIZE],
        veloc_x: &[f32; ARRAY_SIZE],
        veloc_y: &[f32; ARRAY_SIZE],
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

                d[_IX!(i, j, size_n)] = s0
                    * (t0 * d0[_IX!(i0i, j0i, size_n)] + (t1 * (d0[_IX!(i0i, j1i, size_n)])))
                    + s1 * (t0 * d0[_IX!(i1i, j0i, size_n)] + (t1 * (d0[_IX!(i1i, j1i, size_n)])));
                i += 1;
                ifloat += 1.0;
            }
            j += 1;
            jfloat += 1.0;
        }
        set_bnd(b, d, size_n);
    }

    pub fn project(
        veloc_x: &mut [f32; ARRAY_SIZE],
        veloc_y: &mut [f32; ARRAY_SIZE],
        p: &mut [f32; ARRAY_SIZE],
        div: &mut [f32; ARRAY_SIZE],
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
