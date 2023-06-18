#!/bin/bash
bashrc_path="$HOME/.bashrc"
repo_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# make `bin` dir if not exists
mkdir -p bin

# Remove all links in bin
echo "Removing all links in bin..."
find bin -type l -delete


# Copy resources
echo "Copying resources..."
cp -r "resources" "bin/"

# Install bash functions
echo "Installing bash functions..."

bash_src_file="src/bash/cli_utils_bash.sh"

line1="export CLI_UTILS_INSTALL_DIR=\"$repo_dir\""
line2="source $repo_dir/$bash_src_file"
# Add only if not already added
if ! grep -q "$line1" "$bashrc_path"; then
    echo "Adding $line1 to $bashrc_path"
    echo "$line1" >> "$bashrc_path"
fi
if ! grep -q "$line2" "$bashrc_path"; then
    echo "Adding $line2 to $bashrc_path"
    echo "$line2" >> "$bashrc_path"
fi

# Install rust binaries
echo "Installing rust binaries..."

# List the stem of all directories in src/rust/bins, excluding the root directory
all_bins=$(find src/rust/bins -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
cargo build --release
for bin in $all_bins; do
    echo "Installing $bin..."
    ln -sf "../target/release/$bin" "bin/$bin"
done

# Install python scripts
echo "Installing python scripts..."
cp -r src/python/lib bin/
for file in src/python/*.py; do
    echo "Installing $file..."
    base=$(basename "$file")
    ln -sf "../$file" "bin/${base%%.*}" # remove extension
done

# Add bin to PATH
echo "Adding bin to PATH..."
line="export PATH=\"$repo_dir/bin:\$PATH\""
if ! grep -q "$line" "$bashrc_path"; then
    echo "Adding $line to $bashrc_path"
    echo "$line" >> "$bashrc_path"
fi
chmod +x bin/*

echo "Done! Restart your terminal to use the new commands."
