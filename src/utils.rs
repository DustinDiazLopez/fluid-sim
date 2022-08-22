

macro_rules! clamp_meta {
    ($fn_name: ident, $ty:ty) => {
        pub fn $fn_name(v: $ty, min: $ty, max: $ty) -> $ty {
            if v < min {
                min
            } else if v > max {
                max
            } else {
                v
            }
        }
    };
}

clamp_meta!(clamp_f32, f32);
clamp_meta!(clamp_i32, i32);

