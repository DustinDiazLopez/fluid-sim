#!/bin/sh

# TODO: https://bevy-cheatbook.github.io/platforms/wasm.html

rustup target add wasm32-unknown-unknown
rustup target add x86_64-pc-windows-gnu

cargo install trunk
cargo install wasm-bindgen-cli

trunk build

trunk build --release
rm docs/*
mkdir -p docs
cp -r ./dist/* ./docs

#cargo build --release --target x86_64-pc-windows-gnu
