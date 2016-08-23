import pytest
import atomman as am
import numpy as np

class Test_is_int:
    def test_int(self):
        assert am.tools.is_int(1)
        
    def test_long(self):
        assert am.tools.is_int(1L)

    def test_int_list(self):
        assert not am.tools.is_int([1])
               
    def test_np_int32(self):
        assert am.tools.is_int(np.array(1, dtype='int32'))
        
    def test_np_int64(self):
        assert am.tools.is_int(np.array(1, dtype='int64'))
        
    def test_np_int64_list(self):
        assert not am.tools.is_int(np.array([1], dtype='int64'))
        
    def test_bool(self):
        assert not am.tools.is_int(True)
        
    def test_bool_list(self):
        assert not am.tools.is_int([True])

    def test_np_bool(self):
        assert not am.tools.is_int(np.array(False))
        
    def test_np_bool_list(self):
        assert not am.tools.is_int(np.array([False]))
              
    def test_float(self):
        assert not am.tools.is_int(1.0)
        
    def test_float_list(self):
        assert not am.tools.is_int([1.0])

    def test_np_float(self):
        assert not am.tools.is_int(np.array(1.0, dtype='float64'))
        
    def test_np_float_list(self):
        assert not am.tools.is_int(np.array([1.0], dtype='float64'))
        
    def test_str(self):
        assert not am.tools.is_int('1')
        
    def test_np_str(self):
        assert not am.tools.is_int(np.array('1'))     

class Test_is_dtype_int:
    def test_int(self):
        value = np.array(52, dtype=int)
        assert am.tools.is_dtype_int(value.dtype)

    def test_long(self):
        value = np.array(52, dtype=long)
        assert am.tools.is_dtype_int(value.dtype)
        
    def test_int32(self):
        value = np.array(52, dtype='int32')
        assert am.tools.is_dtype_int(value.dtype)

    def test_int64(self):
        value = np.array(52, dtype='int64')
        assert am.tools.is_dtype_int(value.dtype)

    def test_int64_list(self):
        value = np.array([52,124,6], dtype='int64')
        assert am.tools.is_dtype_int(value.dtype)

    def test_bool(self):
        value = np.array(True, dtype=bool)
        assert not am.tools.is_dtype_int(value.dtype)

    def test_bool_list(self):
        value = np.array([True, False], dtype=bool)
        assert not am.tools.is_dtype_int(value.dtype)

    def test_float(self):
        value = np.array(52, dtype='float64')
        assert not am.tools.is_dtype_int(value.dtype)

    def test_str(self):
        value = np.array('52')
        assert not am.tools.is_dtype_int(value.dtype)        
    
        
class Test_is_bool:
    def test_int(self):
        assert not am.tools.is_bool(1)
        
    def test_long(self):
        assert not am.tools.is_bool(1L)

    def test_int_list(self):
        assert not am.tools.is_bool([1])
               
    def test_np_int32(self):
        assert not am.tools.is_bool(np.array(1, dtype='int32'))
        
    def test_np_int64(self):
        assert not am.tools.is_bool(np.array(1, dtype='int64'))
        
    def test_np_int64_list(self):
        assert not am.tools.is_bool(np.array([1], dtype='int64'))
        
    def test_bool(self):
        assert am.tools.is_bool(True)
        
    def test_bool_list(self):
        assert not am.tools.is_bool([True])

    def test_np_bool(self):
        assert am.tools.is_bool(np.array(False))
        
    def test_np_bool_list(self):
        assert not am.tools.is_bool(np.array([False]))
              
    def test_float(self):
        assert not am.tools.is_bool(1.0)
        
    def test_float_list(self):
        assert not am.tools.is_bool([1.0])

    def test_np_float(self):
        assert not am.tools.is_bool(np.array(1.0, dtype='float64'))
        
    def test_np_float_list(self):
        assert not am.tools.is_bool(np.array([1.0], dtype='float64'))
        
    def test_str(self):
        assert not am.tools.is_bool('1')
        
    def test_np_str(self):
        assert not am.tools.is_bool(np.array('1'))     


class Test_is_dtype_bool:
    def test_int(self):
        value = np.array(52, dtype=int)
        assert not am.tools.is_dtype_bool(value.dtype)

    def test_long(self):
        value = np.array(52, dtype=long)
        assert not am.tools.is_dtype_bool(value.dtype)
        
    def test_int32(self):
        value = np.array(52, dtype='int32')
        assert not am.tools.is_dtype_bool(value.dtype)

    def test_int64(self):
        value = np.array(52, dtype='int64')
        assert not am.tools.is_dtype_bool(value.dtype)

    def test_int64_list(self):
        value = np.array([52,124,6], dtype='int64')
        assert not am.tools.is_dtype_bool(value.dtype)

    def test_bool(self):
        value = np.array(True, dtype=bool)
        assert am.tools.is_dtype_bool(value.dtype)

    def test_bool_list(self):
        value = np.array([True, False], dtype=bool)
        assert am.tools.is_dtype_bool(value.dtype)

    def test_float(self):
        value = np.array(52, dtype='float64')
        assert not am.tools.is_dtype_bool(value.dtype)

    def test_str(self):
        value = np.array('52')
        assert not am.tools.is_dtype_bool(value.dtype)                