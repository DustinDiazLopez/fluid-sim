#!/bin/sh

# TODO: https://bevy-cheatbook.github.io/platforms/wasm.html
rustup target add wasm32-unknown-unknown
cargo install wasm-server-runner

#cargo build --release
#cargo build --target wasm32-unknown-unknown --release
#cp ./target/wasm32-unknown-unknown/release/fluid_rs.wasm ./docs/.
