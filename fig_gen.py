from gpx_route_matching import path_statistics
import matplotlib.pyplot as plt
import geopandas as gpd
import contextily as ctx
import gpxpy
import sys
import yaml

from utils import gpx_route_to_gdf
from config import Config
from matplotlib import colormaps

"""
Generate figures regarding route choice
Intended to be ran after gpx_route_matching.py has been run
"""

if __name__ == "__main__":

    """
    sys.argv[1] must be path to valid config yaml file
    """
    # Load in the config file
    with open(sys.argv[1], "r") as f:
        config_data = yaml.safe_load(f)
    config = Config.model_validate(config_data)

    path_taken = config.paths.taken_paths_gpkg
    suggested_path_gpks = [f"{config.paths.alternate_routes_prefix}_{i}.gpkg" for i in range(2)]
    personal_gps_track = config.gpx_tracks.gdf_gpx_path
    gpx_tracks = gpd.read_file(personal_gps_track)
    clipped_swisstlm = config.swisstlm3d.clip_gpkg
    streets_layer = config.swisstlm3d.network_layer
    geopkg = gpd.read_file(clipped_swisstlm, layer=streets_layer)
    hochflueregiongpkg = "temp.gpkg"
    geopkg.to_file(hochflueregiongpkg, driver='GPKG')
    taken_paths = gpd.read_file(path_taken)


    """
    Map GPX track on basemap
    """
    if True:
        fig, ax = plt.subplots()
        gpx_tracks.plot(ax=ax, color="orange")
        ctx.add_basemap(ax, crs=gpx_tracks.crs, source=ctx.providers.SwissFederalGeoportal.NationalMapColor)
        ax.set_title('GPX Route Points Plotted')
        plt.savefig(f"{config.name.name}_gpx_plotted.png")

    """
    Map all the paths with the GPX superimposed
    """
    if True:
        fig, ax = plt.subplots()
        geopkg.plot(ax=ax, color="blue", zorder=1)
        gpx_tracks.plot(ax=ax, color="orange", zorder=2)
        ctx.add_basemap(ax, crs=gpx_tracks.crs, source=ctx.providers.SwissFederalGeoportal.NationalMapColor)
        ax.set_title('GPX path superimposed on path network')
        plt.savefig(f"{config.name.name}_paths_plotted.png")

    """
    Map all the paths with the map-matched path superimposed
    """
    if True:
        fig, ax = plt.subplots()
        geopkg.plot(ax=ax, color="blue", zorder=1)
        taken_paths.plot(ax=ax, color="orange", zorder=2, linewidth=4.5)
        ctx.add_basemap(ax, crs=gpx_tracks.crs, source=ctx.providers.SwissFederalGeoportal.NationalMapColor)
        ax.set_title('Taken route superimposed on path network')
        plt.savefig(f"{config.name.name}_taken_route_over_all_paths.png")

    """
    Generate a figure displaying the taken route and the proposed alternate routes
    by OSRM
    """
    if True:
        # Load in the GPX points and display them
        fig, ax = plt.subplots()
        for a_idx, alt_route in enumerate(suggested_path_gpks):
            route = gpd.read_file(alt_route)
            if a_idx == 0:
                route.plot(ax=ax, color="blue", label="OSRM Routes (3)", linewidth=4, zorder=2)
            else:
                route.plot(ax=ax, color="blue", linewidth=4, zorder=2)
        taken_paths.plot(ax=ax, color="orange", label="Taken Route", linewidth=4, zorder=2)
        geopkg.plot(ax=ax, color="black", linewidth=0.15, zorder=1)
        ctx.add_basemap(ax, crs=taken_paths.crs, source=ctx.providers.SwissFederalGeoportal.NationalMapColor)
        ax.set_title("Taken route vs OSRM alternate routes")
        ax.legend()
        plt.savefig(f"{config.name.photo_loc}_path_vs_alternate_paths.png")
    
    def taken_path_attribute_fig_gen(path_gpkg, attribute, title, fig_filename):
        taken = path_statistics(path_gpkg)
        colors = {label: color  for label, color in zip(sorted(list(taken[attribute].keys())), colormaps['Dark2'].colors)}
        fig, ax = plt.subplots()
        ax.pie(list(taken[attribute].values()),
                labels=tuple(taken[attribute].keys()),
                colors=[colors[v] for v in taken[attribute].keys()])
        ax.set_title(title)
        plt.savefig(fig_filename)
        return colors
    
    def suggested_paths_attribute_fig_gen(suggested_path_gpks, attribute, title, fig_filename, colorway):
        cumulative_characteristics = {}
        for alt_route in suggested_path_gpks:
            stats = path_statistics(alt_route)
            desired_stats = stats[attribute]
            for k in desired_stats:
                if k not in cumulative_characteristics:
                    cumulative_characteristics[k] = desired_stats[k]
                else:
                    cumulative_characteristics[k] += desired_stats[k]
        fig, ax = plt.subplots()
        c_idx = 0
        for t in list(cumulative_characteristics.keys()):
            if t not in colorway:
                colorway[t] = colormaps['Set1'].colors[c_idx % len(colormaps['Set1'].colors)]
                c_idx += 1
        ax.pie(list(cumulative_characteristics.values()),
               labels=tuple(cumulative_characteristics.keys()),
               colors=[colorway[v] for v in cumulative_characteristics.keys()])
        ax.set_title(title)
        plt.savefig(fig_filename)

    """
    Get information about the paths taken
    """
    if True:
        cbelag = taken_path_attribute_fig_gen(path_taken, "belagsart",
                                      f"{config.name.name} trail surface type for taken route",
                                        f"{config.name.photo_loc}_taken_route_surfaces.png")
        cobjektarg = taken_path_attribute_fig_gen(path_taken, "objektart",
                                      f"{config.name.name} type of paths taken", f"{config.name.photo_loc}_taken_path_types.png")
        cwanderwege = taken_path_attribute_fig_gen(path_taken, "wanderwege",
                                      f"{config.name.name} wanderwege Taken", f"{config.name.photo_loc}_taken_wanderwege.png")


    """
    Get infromation about OSRM suggested paths
    """
    if True:
        suggested_paths_attribute_fig_gen(suggested_path_gpks, "belagsart",
                                           f"{config.name.name} trail surface type for OSRM routes",
                                           f"{config.name.photo_loc}_alternate_route_surfaces.png",
                                           cbelag)
        suggested_paths_attribute_fig_gen(suggested_path_gpks, "objektart",
                                           f"{config.name.name} type of paths suggested by OSRM",
                                           f"{config.name.photo_loc}_osrm_path_types.png",
                                           cobjektarg)
        suggested_paths_attribute_fig_gen(suggested_path_gpks, "wanderwege",
                                           f"{config.name.name} wanderwege OSRM", 
                                          f"{config.name.photo_loc}_osrm_wanderwege.png",
                                          cwanderwege)
    
    """
    Get information about ALL the possible paths
    """
    if True:
        suggested_paths_attribute_fig_gen([hochflueregiongpkg], "belagsart",
                                    'Trail surface type for all paths',
                                    f"{config.name.photo_loc}_allpaths_route_surfaces.png",
                                    cbelag)
        suggested_paths_attribute_fig_gen([hochflueregiongpkg], "objektart",
                                    "Type of paths (all paths)",
                                    f"{config.name.photo_loc}_allpaths_path_types.png",
                                    cobjektarg)
        suggested_paths_attribute_fig_gen([hochflueregiongpkg], "wanderwege",
                                     "Wanderwege (All Paths)",
                                     f"{config.name.photo_loc}_allpaths_wanderwege.png",
                                     cwanderwege)
