import json
from dataclasses import dataclass
from enum import Enum
from math import sqrt
from typing import List


ELEMENT_RADII = {
    'C': 1.70,
    'O': 1.52,
    'O1': 1.52,
    'S': 1.80,
    'N': 1.55,
    'N1': 1.55,
    'H': 1.20,
}


class InteractionKind(Enum):
    HYDROPHOBIC = 0
    DONOR = 1
    ACCEPTOR = 2
    EXCLUSION = 3


@dataclass
class Feature:
    type: InteractionKind
    x: float
    y: float
    z: float
    radius: float
    weight: float = 1


def distance(a: Feature, b: Feature) -> float:
    return sqrt((b.x-a.x)**2 + (b.y-a.y)**2 + (b.z-a.z)**2)


def merge_two_features(a: Feature, b: Feature) -> Feature:
    assert a.type == b.type
    return Feature(
        type=a.type,
        x=(a.x*a.weight + b.x*b.weight)/(a.weight+b.weight),
        y=(a.y*a.weight + b.y*b.weight)/(a.weight+b.weight),
        z=(a.z*a.weight + b.z*b.weight)/(a.weight+b.weight),
        radius=(a.radius*a.weight + b.radius*b.weight)/(a.weight+b.weight),
        weight=a.weight+b.weight,
    )


def maybe_merge_nearby_features(features: [Feature], threshold: float = 2) -> List[Feature]:
    """Mescla features muito próximos pela média ponderada.

    """
    all_features = features[:]
    dirty_pass = False
    while not dirty_pass and len(all_features) >= 2:
        new_features = []
        for i, feat1 in enumerate(all_features[:]):
            was_merged = False
            for j, feat2 in enumerate(all_features[:]):
                if i == j:
                    dirty_pass = True
                    continue
                if feat1.type == feat2.type and \
                        distance(feat1, feat2) <= threshold:
                    all_features.remove(feat1)
                    all_features.remove(feat2)
                    new_features.append(merge_two_features(feat1, feat2))
                    was_merged = True
                    dirty_pass = False
                    break
                else:
                    dirty_pass = False
            if was_merged:
                break
        all_features.extend(new_features)
    return all_features


class PharmacophoreJsonWriter:
    """Write a list of features to a Pharmit/ZINCPharmer session file.

    """

    def _build_feature(self, feat: Feature):
        if feat.type == InteractionKind.HYDROPHOBIC:
            return {
                "name": "Hydrophobic",
                "enabled": True,
                "x": feat.x,
                "y": feat.y,
                "z": feat.z,
                "radius": feat.radius,
                "weight": feat.weight,
            }
        elif feat.type == InteractionKind.DONOR:
            return {
                "name": "HydrogenDonor",
                "enabled": True,
                "x": feat.x,
                "y": feat.y,
                "z": feat.z,
                "radius": feat.radius,
                "weight": feat.weight,
            }
        elif feat.type == InteractionKind.ACCEPTOR:
            return {
                "name": "HydrogenAcceptor",
                "enabled": True,
                "x": feat.x,
                "y": feat.y,
                "z": feat.z,
                "radius": feat.radius,
                "weight": feat.weight,
            }
        elif feat.type == InteractionKind.EXCLUSION:
            return {
                "name": "ExclusionSphere",
                "enabled": True,
                "x": feat.x,
                "y": feat.y,
                "z": feat.z,
                "radius": feat.radius,
                "weight": feat.weight,
            }
        else:
            raise Exception("Unexpected error")

    def write(self, features: [Feature], filename: str) -> None:
        with open(filename, 'w') as fp:
            return json.dump({
                'points': [
                    self._build_feature(feat) for feat in features
                ]
            }, fp)


@dataclass
class Cluster:
    type: InteractionKind
    atoms: [(float, float, float, str)]
    cluster_size: int


class FakeLigandReader:

    def read(self, filename_rings: str, filename_donors: str, filename_acceptors: str) -> [Feature]:
        clusters = []
        clusters.extend(self._read_sdf(
            filename_rings, InteractionKind.HYDROPHOBIC))
        clusters.extend(self._read_sdf(filename_donors, InteractionKind.DONOR))
        clusters.extend(self._read_sdf(
            filename_acceptors, InteractionKind.ACCEPTOR))

        features = []
        for cluster in clusters:
            features.append(self._cluster_to_feature(cluster))
        return features


    def _read_sdf(self, filename: str, cluster_type: InteractionKind) -> List[Cluster]:
        clusters = []
        with open(filename) as file:
            lines = file.readlines()
            atoms = []
            cluster_size = None

            while len(lines) > 0:
                line = lines.pop(0)
                if line.startswith('   '):
                    chunks = line.split()
                    x = float(chunks[0])
                    y = float(chunks[1])
                    z = float(chunks[2])
                    elem = chunks[3]
                    xyze = (x, y, z, elem)
                    atoms.append(xyze)
                elif line.startswith('> <CLUSTER_SIZE>'):
                    cluster_size = int(lines.pop(0))
                elif line.startswith('$$$$'):
                    cluster = Cluster(
                        type=cluster_type,
                        atoms=atoms,
                        cluster_size=cluster_size
                    )
                    clusters.append(cluster)
                    atoms = []
        return clusters

    def _cluster_to_feature(self, cluster: Cluster) -> Feature:
        if cluster.type == InteractionKind.HYDROPHOBIC:
            x_tot = 0
            y_tot = 0
            z_tot = 0
            for atom in cluster.atoms:
                x, y, z, elem = atom
                x_tot += x
                y_tot += y
                z_tot += z
            return Feature(
                type=InteractionKind.HYDROPHOBIC,
                x=x_tot / len(cluster.atoms),
                y=y_tot / len(cluster.atoms),
                z=z_tot / len(cluster.atoms),
                radius=1,
                weight=cluster.cluster_size,
            )
        elif cluster.type == InteractionKind.DONOR:
            assert len(cluster.atoms) == 2
            assert cluster.atoms[0][3] == 'N'
            assert cluster.atoms[1][3] == 'H'
            x, y, z, elem = cluster.atoms[0]
            return Feature(
                type=InteractionKind.DONOR,
                x=x,
                y=y,
                z=z,
                radius=1,
                weight=cluster.cluster_size
            )
        elif cluster.type == InteractionKind.ACCEPTOR:
            assert len(cluster.atoms) == 1
            assert cluster.atoms[0][3] == 'O'
            x, y, z, elem = cluster.atoms[0]
            return Feature(
                type=InteractionKind.ACCEPTOR,
                x=x,
                y=y,
                z=z,
                radius=1,
                weight=cluster.cluster_size
            )



class Strategy:

    def run(feat: List[Feature]) -> List[Feature]:
        raise NotImplementedError


class SimpleStrategy(Strategy):

    def run(feats: List[Feature], num_points: int = 5) -> List[Feature]:
        return maybe_merge_nearby_features(feats)[:num_points]


class Adjusted(Strategy):

    def run(feats: List[Feature],
            num_points: int = 5,
            min_radius: float = 0,
            max_radius: float = 1
    ) -> List[Feature]:
        pass


def find_exclusion_features(filename_site_pdb: str) -> List[Feature]:
    feats: List[Feature] = []
    with open(filename_site_pdb) as pdb_fp:
        for line in pdb_fp:
            if line.startswith('ATOM'):
                elem = line[77:79]
                radius = ELEMENT_RADII[elem.strip()]
                feat = Feature(
                    InteractionKind.EXCLUSION,
                    x=float(line[31:39]),
                    y=float(line[39:47]),
                    z=float(line[47:55]),
                    radius=radius
                )
                feats.append(feat)
    return feats


def fake2json(
    filename_rings: str,
    filename_donors: str,
    filename_acceptors: str,
    filename_site_pdb: str,
    filename_output: str,
) -> None:


    reader = FakeLigandReader()
    feats = reader.read(filename_rings, filename_donors, filename_acceptors)
    new_feats = maybe_merge_nearby_features(feats)
    writer = PharmacophoreJsonWriter()
    site = find_exclusion_features(filename_site_pdb)

    top_feats = list(sorted(new_feats, key=lambda f: -f.weight))
    from pprint import pprint
    all_feats = top_feats + site
    pprint(all_feats)
    writer.write(all_feats, filename_output)

# TODO distância máxima (8A) entre pontos
# TODO filtrar pontos pelo DC
