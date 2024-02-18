import os
import tempfile
import unittest

import numpy as np
import tensorflow

from medaka import models
from medaka.common import Sample
from medaka.datastore import DataStore
from medaka.features import BaseFeatureEncoder
from medaka.labels import BaseLabelScheme
import medaka.options


class TestModelFiles(unittest.TestCase):

    def test_000_total_bundled_size(self):
        total = 0
        for name in medaka.options.current_models:
            model_file = models.resolve_model(name)
            total += os.path.getsize(model_file)
        self.assertLess(total / 1024 / 1024, 45, "Bundled model file size too large")

    def test_001_default_models(self):
        for name in medaka.options.default_models.values():
            if name not in medaka.options.current_models:
                self.fail('Default Model {} not in current_models'.format(model_file))

    def test_010_failed_download(self):
        name = 'garbage'
        medaka.options.allowed_models.append(name)
        with self.assertRaises(medaka.models.DownloadError):
            models.resolve_model(name)
        medaka.options.allowed_models.pop()

    def test_011_success_download(self):
        name = 'r941_min_high_g351'
        model_file = models.resolve_model(name)
        tmp_file = "{}.tmp".format(model_file)
        os.rename(model_file, tmp_file)
        new_file = models.resolve_model(name)
        self.assertTrue(os.path.isfile(new_file))
        os.remove(new_file)
        os.rename(tmp_file, model_file)

    def test_999_load_all_models(self):
        for name in medaka.options.allowed_models:
            model_file = models.resolve_model(name)
            model = medaka.models.open_model(model_file).load_model()
            self.assertIsInstance(model, tensorflow.keras.models.Model)
            # Check we can get necessary functions for inference
            with DataStore(model_file) as ds:
                feature_encoder = ds.get_meta('feature_encoder')
                self.assertIsInstance(feature_encoder, BaseFeatureEncoder)
                label_scheme = ds.get_meta('label_scheme')
                self.assertIsInstance(label_scheme, BaseLabelScheme)


class TestBuildModel(unittest.TestCase):

    def test_000_build_all_models(self):
        num_classes, time_steps, feat_len = 5, 5, 5
        for name, func in models.model_builders.items():
            model = func(feat_len, num_classes, time_steps=time_steps)



class TestMajorityModel(unittest.TestCase):

    def test_000_initialise_majority_model(self):
        majority_model = models.build_majority(10, 5)
        self.assertIsInstance(majority_model, tensorflow.keras.models.Model)
