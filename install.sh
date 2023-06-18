#!/bin/bash

vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}

bashrc_path="$HOME/.bashrc"
repo_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# If rustc is not installed, error
if ! command -v rustc &> /dev/null
then
    echo "rustc could not be found. Please install rustc."
    exit
fi

# Check rustc version
echo "Checking rustc version..."
rustc_version=$(rustc --version | awk '{print $2}')
# if version older than 1.64.0, update
vercomp "$rustc_version" "1.64.0"
if [ $? = 2 ]
then
    echo "rustc version is $rustc_version. Please update rustc to 1.64.0 or newer."
    exit
fi


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
