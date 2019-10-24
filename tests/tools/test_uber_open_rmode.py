# coding: utf-8

# Standard Python libraries
import os

# https://docs.pytest.org/en/latest/
import pytest

from atomman.tools import uber_open_rmode

class Test_uber_open_rmode(object):

    @property
    def content(self):
        return "This is the content of my file."

    def test_stringread(self):
        with uber_open_rmode(self.content) as f:
            content = f.read()
        assert content.decode('UTF-8') == self.content

    def test_fileread(self, tmpdir):
        contentfile = os.path.join(str(tmpdir), 'content.txt')
        with open(contentfile, 'w') as f:
            f.write(self.content)
        
        with uber_open_rmode(contentfile) as f:
            content = f.read()
        assert content.decode('UTF-8') == self.content

    def test_objectread(self, tmpdir):
        contentfile = os.path.join(str(tmpdir), 'content.txt')
        with open(contentfile, 'w') as f:
            f.write(self.content)

        with open(contentfile, 'rb') as openf:
            with uber_open_rmode(openf) as f:
                content = f.read()
        assert content.decode('UTF-8') == self.content
