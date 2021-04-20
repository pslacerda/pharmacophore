from pharmacophore import FakeLigandReader, merge_two_features, Feature, InteractionKind, maybe_merge_nearby_features


def test_merge_two_features():
    result = merge_two_features(
        Feature(x=1, y=1, z=1, radius=1, weight=1, type=InteractionKind.DONOR),
        Feature(x=2, y=2, z=2, radius=1, weight=1, type=InteractionKind.DONOR)
    )
    expected = Feature(x=1.5, y=1.5, z=1.5, radius=1,
                       weight=2, type=InteractionKind.DONOR)
    assert result == expected


def test_maybe_merge_nearby_features():
    result = maybe_merge_nearby_features([
        Feature(x=1, y=1, z=1, radius=1, weight=1, type=InteractionKind.DONOR),
        Feature(x=2, y=2, z=2, radius=1, weight=1, type=InteractionKind.DONOR),
        Feature(x=1.5, y=1.5, z=1.5, radius=1, weight=1, type=InteractionKind.DONOR),
        Feature(x=9, y=9, z=9, radius=1, weight=1, type=InteractionKind.DONOR),
    ])
    expected = [
        Feature(x=1.5, y=1.5, z=1.5, radius=1,
                weight=3, type=InteractionKind.DONOR),
        Feature(x=9, y=9, z=9, radius=1, weight=1, type=InteractionKind.DONOR)
    ]
    equals = 0
    for feat1 in result:
        for feat2 in result:
            if feat1 == feat2:
                equals += 1
    assert equals == 2


def test_FakeLigandReader():
    r = FakeLigandReader()
    clusters = r._read_sdf(
        'examples/3cqw/v2_3cqw_cutoff_50_donors_medchem.sdf', InteractionKind.DONOR)
    assert len(clusters) == 8

    feat = r._cluster_to_feature(clusters[0])
    assert feat.x == 40.4813

    feat_list = r.read('examples/3cqw/v2_3cqw_cutoff_50_rings.sdf',
                       'examples/3cqw/v2_3cqw_cutoff_50_donors_medchem.sdf',
                       'examples/3cqw/v2_3cqw_cutoff_50_acceptors_medchem.sdf')
    assert isinstance(feat_list[0], Feature)


if __name__ == '__main__':
    test_merge_two_features()
    test_maybe_merge_nearby_features()
    test_FakeLigandReader()
    print('OK')
