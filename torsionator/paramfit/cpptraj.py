import os


def write_xyz_to_nc_in(scan_dir: str, leap_dir: str, method_subdir: str) -> None:
    """Write a cpptraj input file to convert geometries.xyz to geometries.nc."""
    out_dir = os.path.join(scan_dir, method_subdir)
    parm_path = os.path.join(scan_dir, leap_dir, "MOL.top")
    trajin = os.path.join(out_dir, "geometries.xyz")

    with open(os.path.join(out_dir, "xyz_to_nc.in"), "w") as f:
        f.write(f"parm {parm_path}\n")
        f.write(f"trajin {trajin}\n")
        f.write(f"parmbox parm {parm_path} x 100 y 100 z 100 alpha 90 beta 90 gamma 90\n")
        f.write("box x 100 y 100 z 100 alpha 90 beta 90 gamma 90\n")
        f.write("trajout geometries.nc netcdf\n")
