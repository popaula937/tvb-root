import os
import re

import numpy
import siibra


### Extra dependencies: siibra

def example_notebook():
    atlas = siibra.atlases["human"]

    # in MNI 152 space
    atlas.select(parcellation="julich")
    icbm_map = atlas.get_map(space="mni152")
    # the julich brain map comes in separate l/r hemispheres,
    # so we iterate over all maps.
    # for m in icbm_map.fetchall():
    # plotting.plot_stat_map(m)

    # bigbrain
    reso_mm = 0.64
    bigbrain_tpl = atlas.get_template("bigbrain")
    bigbrain_map = atlas.get_map(space="bigbrain")
    # plotting.plot_stat_map(bigbrain_map.fetch(reso_mm),
    #                        bigbrain_tpl.fetch(reso_mm))

    # DK atlas
    atlas.select(parcellation="desikan")
    dk_map = atlas.get_map(space="mni152")
    # plotting.plot_stat_map(dk_map.fetch(), cmap=plt.cm.tab10)


def get_connectivity_matrices(connectivities):
    for conn in connectivities:
        # Averaged_SC_JuBrain_246Regions_eNKI_10M_count_MEAN
        # Averaged_SC_JuBrain_246Regions_HCP_10M_count_MEAN
        if 'Averaged_SC_JuBrain_246Regions_HCP_10M_count_MEAN' in conn.name:
            conn_weights = conn
            print("Weights matrix: {}".format(conn_weights.name))

        # Averaged_SC_JuBrain_246Regions_eNKI_10M_length_MEAN
        # Averaged_SC_JuBrain_246Regions_HCP_10M_length_MEAN
        if 'Averaged_SC_JuBrain_246Regions_HCP_10M_length_MEAN' in conn.name:
            conn_tracts = conn
            print("Tracts matrix: {}".format(conn_tracts.name))

    return conn_weights, conn_tracts


def get_position(position_list, idx):
    region_position = position_list[idx].get('position')
    if region_position is None:
        region_position = [0, 0, 0]
    else:
        region_position = [reg_pos / 1000000 for reg_pos in region_position]
    return region_position


def example_vep():
    atlas = siibra.atlases["human"]
    parcellation = siibra.parcellations.JULICH_BRAIN_CYTOARCHITECTONIC_MAPS_2_5
    jubrain = atlas.get_parcellation(parcellation)
    # julichbrain_mpm_left = atlas.get_map('mni152', maptype="labelled").fetch()
    # targetspace = atlas.get_space('mni152')

    mpm = atlas.get_map(space="mni152", parcellation="julich")

    # If I use the jubrain parcellation 2.9 for querying the connectivities, I get a single one as result.
    # If I use the atlas for querying, I get multiple connectivities back, but the logs say that the 2.9 julich parcellation was used by default?
    connectivities = siibra.get_features(jubrain, siibra.modalities.ConnectivityMatrix)

    conn_weights, conn_tracts = get_connectivity_matrices(connectivities)
    reg_labels = conn_weights.regionnames
    weights = conn_weights.array
    tracts = conn_tracts.array

    regex = '([\w\s.-]+)( \(([\w\s,]+)\))? (left|right|both)?'
    p = re.compile(regex)
    # p.match('Area 4p (PreCG) - left hemisphere')

    reg_labels_corrected = list()
    reg_names_corrected = list()
    hemi_array = list()
    reg_positions = list()
    for reg_label in reg_labels:
        match_rez = p.match(reg_label)
        if match_rez:
            region_name = match_rez.groups()[0]
            region_label = match_rez.groups()[-2]
            region_hemi = match_rez.groups()[-1]

            if region_hemi == 'both':
                continue

            reg_labels_corrected.append(region_label)
            reg_names_corrected.append(region_name.replace(' ', '_'))

            region = atlas.get_region(region_name, parcellation)
            region_position_list = region.attrs.get('children')

            if region_hemi == 'right':
                hemi_array.append(1)
                if len(region_position_list) > 1:
                    region_position = get_position(region_position_list, 1)
                else:
                    region_position = get_position(region_position_list, 0)
            else:
                hemi_array.append(0)
                region_position = get_position(region_position_list, 0)

            reg_positions.append(region_position)

    centres_content = numpy.concatenate((numpy.array(reg_names_corrected)[:, None], reg_positions), axis=1)

    out_dir = 'julich_conn_pos_altered'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    numpy.savetxt(os.path.join(out_dir, "julich_centers.txt"), centres_content, "%s")
    numpy.savetxt(os.path.join(out_dir, "julich_hemispheres.txt"), hemi_array, "%s")
    numpy.savetxt(os.path.join(out_dir, "julich_weights.txt"), weights, "%f")
    numpy.savetxt(os.path.join(out_dir, "julich_tracts.txt"), tracts, "%f")

    # atlas.clear_selection()
    # features_parc = atlas.get_features(siibra.modalities.GeneExpression, gene=siibra.features.gene_names.GABARAPL2)
    # features_parc = atlas.get_features(siibra.modalities.ConnectivityMatrix)
    # print(features_parc)


os.environ[
    'HBP_AUTH_TOKEN'] = 'eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICJfNkZVSHFaSDNIRmVhS0pEZDhXcUx6LWFlZ3kzYXFodVNJZ1RXaTA1U2k0In0.eyJleHAiOjE2MzE2MTkzMzUsImlhdCI6MTYzMTAxNDUzNSwiYXV0aF90aW1lIjoxNjMxMDE0NTM1LCJqdGkiOiJmNjI5NTc1OS04ZDU3LTQzYmMtYTgwMC0yYTgwYTlhMTMwYmQiLCJpc3MiOiJodHRwczovL2lhbS5lYnJhaW5zLmV1L2F1dGgvcmVhbG1zL2hicCIsImF1ZCI6WyJyZWFsbS1tYW5hZ2VtZW50IiwieHdpa2kiLCJ0ZWFtIiwiZ3JvdXAiXSwic3ViIjoiMzY3MjQ4OWEtOGI4ZS00NGFkLWI4M2ItNTQ2MzNjYjVhNzVlIiwidHlwIjoiQmVhcmVyIiwiYXpwIjoianVweXRlcmh1YiIsInNlc3Npb25fc3RhdGUiOiJlYjY1YWFhMS0yOGYwLTQzM2ItOTAyYy1kYzY5NmY5MjYyMGYiLCJhY3IiOiIxIiwiYWxsb3dlZC1vcmlnaW5zIjpbImh0dHBzOi8vanVweXRlcmh1Yi5hcHBzLmpzYy5oYnAuZXUvIiwiaHR0cHM6Ly9sYWIuZWJyYWlucy5ldS8iLCJodHRwczovL2xhYi5qc2MuZWJyYWlucy5ldS8iXSwicmVhbG1fYWNjZXNzIjp7InJvbGVzIjpbIm9mZmxpbmVfYWNjZXNzIl19LCJzY29wZSI6InByb2ZpbGUgY29sbGFiLmRyaXZlIG9mZmxpbmVfYWNjZXNzIGNsYi53aWtpLndyaXRlIGNsYi53aWtpLnJlYWQgdGVhbSBjbGIuZHJpdmU6d3JpdGUgcm9sZXMgZW1haWwgb3BlbmlkIGNsYi5kcml2ZTpyZWFkIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsIm5hbWUiOiJQYXVsYSBQb3BhIiwibWl0cmVpZC1zdWIiOiIzMDc0NDEiLCJwcmVmZXJyZWRfdXNlcm5hbWUiOiJwYXVsYXBvcGEiLCJnaXZlbl9uYW1lIjoiUGF1bGEiLCJmYW1pbHlfbmFtZSI6IlBvcGEiLCJlbWFpbCI6InBhdWxhLnBvcGFAY29kZW1hcnQucm8ifQ.p6AsCUgSDJmVVqL2adNNW9BGEZPmMc-r9socuIFVLorfQllx2-cpy3v3BCWw2STTeKvfreqbaK1uRvc1fJpm441jNBCo5Hhy8UohMp_9HrwiyDHU4YXttmTc8fCkwhANwfDHA3VuVukHJUYpKjOn7_MsYZjjKE04MM94z9_iOplT9VzYCtG3gMnNfzPYkuelLAudcgojSGqjvcdu74H0nrDlYSGLfAEADtjpkdOY8V0UASioFYmduxAn8Zp_T5_wqcipradxbp1z2vpriFJzQuVPahw_k32B-AHdEgnFhfBJQSZuO45J62zK2xTU1kVkA90DDXnC5CIwObFLg8pE1Q'
example_vep()

# "Averaged_SC_JuBrain_184Regions_HCP_10M_count_MEAN"
# "Averaged_SC_JuBrain_184Regions_HCP_10M_length_MEAN"
