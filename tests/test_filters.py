#!/bin/env python

import unittest
import filters


class TestFilters(unittest.TestCase):

    def setUp(self):
        self.text_a = "key_a:<1.0,2.0>"
        self.text_b = "key_b:<1,2>"
        self.text_c = "key_c:<1,->"
        self.text_d = "key_d:<-,2>"
        self.text_e = "key_e:SomeString"
        self.text_f = "key_f:[a,b,c]"
        self.text_g = "key_g~String"
        self.text_h = "key_h!:SomeString"
        self.text_i = "key_i!~SomeString"

    def runTest(self):
        factory = filters.FilterFactory()

        # numeric a
        key, flt = factory.get_filter(self.text_a)
        self.assertEqual(key, "key_a")
        self.assertEqual(flt._minimum, 1.0)
        self.assertEqual(flt._maximum, 2.0)
        with self.assertRaises(filters.RejectedException):
            flt.check(3)
        with self.assertRaises(filters.RejectedException):
            flt.check(0.5)

        x = False
        try:
            flt.check(1.5)
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # numeric b
        key, flt = factory.get_filter(self.text_b)
        self.assertEqual(key, "key_b")
        self.assertEqual(flt._minimum, 1)
        self.assertEqual(flt._maximum, 2)
        with self.assertRaises(filters.RejectedException):
            flt.check(3)
        with self.assertRaises(filters.RejectedException):
            flt.check(0.5)

        x = False
        try:
            flt.check(1.5)
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # numeric c
        key, flt = factory.get_filter(self.text_c)
        self.assertEqual(key, "key_c")
        self.assertEqual(flt._minimum, 1)
        self.assertEqual(flt._maximum, None)
        with self.assertRaises(filters.RejectedException):
            flt.check(0.5)

        x = False
        try:
            flt.check(1.5)
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # numeric d
        key, flt = factory.get_filter(self.text_d)
        self.assertEqual(key, "key_d")
        self.assertEqual(flt._minimum, None)
        self.assertEqual(flt._maximum, 2)
        with self.assertRaises(filters.RejectedException):
            flt.check(3)

        x = False
        try:
            flt.check(1.5)
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # string e
        key, flt = factory.get_filter(self.text_e)
        self.assertEqual(key, "key_e")
        self.assertEqual(flt._equals, "SomeString")
        with self.assertRaises(filters.RejectedException):
            flt.check("String")

        x = False
        try:
            flt.check("SomeString")
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # set f
        key, flt = factory.get_filter(self.text_f)
        self.assertEqual(key, "key_f")
        self.assertEqual(flt._in, ["a", "b", "c"])
        with self.assertRaises(filters.RejectedException):
            flt.check("d")

        x = False
        try:
            flt.check("b")
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # regexp g
        key, flt = factory.get_filter(self.text_g)
        self.assertEqual(key, "key_g")
        with self.assertRaises(filters.RejectedException):
            flt.check("d")

        x = False
        try:
            flt.check("The word String is contained in this sentence")
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # regexp h
        key, flt = factory.get_filter(self.text_h)
        self.assertEqual(key, "key_h")
        with self.assertRaises(filters.RejectedException):
            flt.check("SomeString")

        # check that SomeString is not Something
        x = False
        try:
            flt.check("Something")
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, True)

        # regexp i
        key, flt = factory.get_filter(self.text_i)
        self.assertEqual(key, "key_i")
        with self.assertRaises(filters.RejectedException):
            flt.check("SomeString")

        x = False
        try:
            flt.check("The word SomeString is contained in this sentence")
            x = True
        except filters.RejectedException:
            pass
        self.assertEqual(x, False)

