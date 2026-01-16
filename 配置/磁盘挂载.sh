sudo fuser -kvm /mnt/windowsdata
sudo umount -f /mnt/windowsdata
UUID=86A233BAA233AE15 /mnt/windowsdata ntfs-3g defaults,uid=1001,gid=1001,dmask=000,fmask=000 0 0
sudo ntfsfix /dev/nvme0n1p3
sudo mount -a
ls -ld /mnt/windowsdata
