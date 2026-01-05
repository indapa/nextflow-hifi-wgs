#!/bin/bash

# Create directories
echo "Creating directories..."
sudo mkdir -p /mnt/work
sudo mkdir -p /mnt/scratch

# Set permissions (Read/Write/Execute for everyone)
echo "Setting permissions..."
sudo chmod -R 777 /mnt/work /mnt/scratch

echo "System setup complete."
