from __future__ import annotations

import ctypes
import shutil
import subprocess
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
BUILD_DIR = ROOT / "build"
SOURCE = ROOT / "cuda_cones.cu"
PTX = BUILD_DIR / "cuda_cones.ptx"


def _needs_rebuild() -> bool:
    if not PTX.exists():
        return True
    try:
        return SOURCE.stat().st_mtime > PTX.stat().st_mtime
    except FileNotFoundError:
        return True


def build_cuda_ptx() -> Path:
    nvcc = shutil.which("nvcc")
    if nvcc is None:
        raise RuntimeError("nvcc not found")
    if not SOURCE.exists():
        raise RuntimeError(f"missing CUDA source: {SOURCE}")

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        nvcc,
        "-ptx",
        "-O3",
        "-arch=compute_52",
        str(SOURCE),
        "-o",
        str(PTX),
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)
    return PTX


def _load_ptx() -> bytes:
    if _needs_rebuild():
        build_cuda_ptx()
    return PTX.read_bytes()


class CUDAError(RuntimeError):
    pass


def _load_driver() -> ctypes.CDLL:
    for name in ("/usr/lib/wsl/lib/libcuda.so.1", "libcuda.so.1", "libcuda.so"):
        try:
            return ctypes.CDLL(name)
        except OSError:
            continue
    raise RuntimeError("libcuda.so not found")


class CUDAEnergyBackend:
    def __init__(self, n: int, dim: int) -> None:
        if n <= 0 or dim <= 0:
            raise ValueError("n and dim must be positive")
        self.n = n
        self.dim = dim
        self.lib = _load_driver()
        self._configure_api()

        self._check(self.lib.cuInit(0), "cuInit")
        self.device = ctypes.c_int()
        self._check(self.lib.cuDeviceGet(ctypes.byref(self.device), 0), "cuDeviceGet")
        self.context = ctypes.c_void_p()
        self._check(self._cu_ctx_create(ctypes.byref(self.context), 0, self.device), "cuCtxCreate")

        self.module = ctypes.c_void_p()
        ptx = _load_ptx()
        ptx_buf = ctypes.create_string_buffer(ptx + b"\0")
        self._check(self.lib.cuModuleLoadData(ctypes.byref(self.module), ctypes.cast(ptx_buf, ctypes.c_void_p)), "cuModuleLoadData")

        self.function = ctypes.c_void_p()
        self._check(
            self.lib.cuModuleGetFunction(ctypes.byref(self.function), self.module, b"compute_energy_kernel"),
            "cuModuleGetFunction",
        )

        self.d_z = ctypes.c_uint64()
        self.d_x = ctypes.c_uint64()
        self.d_r = ctypes.c_uint64()
        self.d_out = ctypes.c_uint64()

        self._check(self._cu_mem_alloc(ctypes.byref(self.d_z), n * n), "cuMemAlloc(d_z)")
        self._check(self._cu_mem_alloc(ctypes.byref(self.d_x), n * dim * ctypes.sizeof(ctypes.c_double)), "cuMemAlloc(d_x)")
        self._check(self._cu_mem_alloc(ctypes.byref(self.d_r), n * ctypes.sizeof(ctypes.c_double)), "cuMemAlloc(d_r)")
        self._check(self._cu_mem_alloc(ctypes.byref(self.d_out), n * n * ctypes.sizeof(ctypes.c_double)), "cuMemAlloc(d_out)")

        self._z_buffer = (ctypes.c_uint8 * (n * n))()
        self._x_buffer = (ctypes.c_double * (n * dim))()
        self._r_buffer = (ctypes.c_double * n)()
        self._out_buffer = (ctypes.c_double * (n * n))()

    def _configure_api(self) -> None:
        self.lib.cuInit.argtypes = [ctypes.c_uint]
        self.lib.cuInit.restype = ctypes.c_int
        self.lib.cuDeviceGet.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        self.lib.cuDeviceGet.restype = ctypes.c_int
        self.lib.cuDeviceGetCount.argtypes = [ctypes.POINTER(ctypes.c_int)]
        self.lib.cuDeviceGetCount.restype = ctypes.c_int
        self.lib.cuModuleLoadData.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p]
        self.lib.cuModuleLoadData.restype = ctypes.c_int
        self.lib.cuModuleGetFunction.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_char_p]
        self.lib.cuModuleGetFunction.restype = ctypes.c_int
        self.lib.cuMemAlloc_v2.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t]
        self.lib.cuMemAlloc_v2.restype = ctypes.c_int
        self.lib.cuMemcpyHtoD_v2.argtypes = [ctypes.c_uint64, ctypes.c_void_p, ctypes.c_size_t]
        self.lib.cuMemcpyHtoD_v2.restype = ctypes.c_int
        self.lib.cuMemcpyDtoH_v2.argtypes = [ctypes.c_void_p, ctypes.c_uint64, ctypes.c_size_t]
        self.lib.cuMemcpyDtoH_v2.restype = ctypes.c_int
        self.lib.cuLaunchKernel.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_void_p),
            ctypes.POINTER(ctypes.c_void_p),
        ]
        self.lib.cuLaunchKernel.restype = ctypes.c_int
        self.lib.cuCtxSynchronize.argtypes = []
        self.lib.cuCtxSynchronize.restype = ctypes.c_int
        self.lib.cuCtxDestroy_v2.argtypes = [ctypes.c_void_p]
        self.lib.cuCtxDestroy_v2.restype = ctypes.c_int
        self.lib.cuModuleUnload.argtypes = [ctypes.c_void_p]
        self.lib.cuModuleUnload.restype = ctypes.c_int
        self.lib.cuGetErrorString.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        self.lib.cuGetErrorString.restype = ctypes.c_int
        self._cu_ctx_create = self.lib.cuCtxCreate_v2 if hasattr(self.lib, "cuCtxCreate_v2") else self.lib.cuCtxCreate
        self._cu_ctx_create.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_uint, ctypes.c_int]
        self._cu_ctx_create.restype = ctypes.c_int
        self._cu_mem_alloc = self.lib.cuMemAlloc_v2 if hasattr(self.lib, "cuMemAlloc_v2") else self.lib.cuMemAlloc
        self._cu_mem_alloc.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_size_t]
        self._cu_mem_alloc.restype = ctypes.c_int

    def _check(self, code: int, where: str) -> None:
        if code == 0:
            return
        msg = ctypes.c_char_p()
        if hasattr(self.lib, "cuGetErrorString") and self.lib.cuGetErrorString(code, ctypes.byref(msg)) == 0 and msg.value:
            raise CUDAError(f"{where}: {msg.value.decode('utf-8', 'replace')}")
        raise CUDAError(f"{where}: CUDA error {code}")

    def set_z(self, z: Sequence[Sequence[bool]]) -> None:
        idx = 0
        for row in z:
            for cell in row:
                self._z_buffer[idx] = 1 if cell else 0
                idx += 1
        self._check(
            self.lib.cuMemcpyHtoD_v2(self.d_z, ctypes.cast(self._z_buffer, ctypes.c_void_p), self.n * self.n),
            "cuMemcpyHtoD(d_z)",
        )

    def compute(self, x: Sequence[Sequence[float]], r: Sequence[float], rave: float) -> list[float]:
        idx = 0
        for row in x:
            for value in row:
                self._x_buffer[idx] = float(value)
                idx += 1
        for i, value in enumerate(r):
            self._r_buffer[i] = float(value)

        self._check(
            self.lib.cuMemcpyHtoD_v2(self.d_x, ctypes.cast(self._x_buffer, ctypes.c_void_p), self.n * self.dim * ctypes.sizeof(ctypes.c_double)),
            "cuMemcpyHtoD(d_x)",
        )
        self._check(
            self.lib.cuMemcpyHtoD_v2(self.d_r, ctypes.cast(self._r_buffer, ctypes.c_void_p), self.n * ctypes.sizeof(ctypes.c_double)),
            "cuMemcpyHtoD(d_r)",
        )

        n_c = ctypes.c_int(self.n)
        dim_c = ctypes.c_int(self.dim)
        rave_c = ctypes.c_double(rave)
        args = (ctypes.c_void_p * 7)(
            ctypes.cast(ctypes.byref(self.d_z), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(self.d_x), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(self.d_r), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(rave_c), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(self.d_out), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(n_c), ctypes.c_void_p),
            ctypes.cast(ctypes.byref(dim_c), ctypes.c_void_p),
        )
        total = self.n * self.n
        threads = 256
        blocks = (total + threads - 1) // threads
        self._check(
            self.lib.cuLaunchKernel(
                self.function,
                blocks,
                1,
                1,
                threads,
                1,
                1,
                0,
                None,
                args,
                None,
            ),
            "cuLaunchKernel",
        )
        self._check(self.lib.cuCtxSynchronize(), "cuCtxSynchronize")
        self._check(
            self.lib.cuMemcpyDtoH_v2(ctypes.cast(self._out_buffer, ctypes.c_void_p), self.d_out, total * ctypes.sizeof(ctypes.c_double)),
            "cuMemcpyDtoH(d_out)",
        )
        return [float(self._out_buffer[i]) for i in range(total)]

    def close(self) -> None:
        if getattr(self, "module", None):
            try:
                self.lib.cuModuleUnload(self.module)
            except Exception:
                pass
            self.module = None
        if getattr(self, "context", None):
            try:
                self.lib.cuCtxDestroy_v2(self.context)
            except Exception:
                pass
            self.context = None

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass


def cuda_device_count() -> int:
    lib = _load_driver()
    lib.cuInit.argtypes = [ctypes.c_uint]
    lib.cuInit.restype = ctypes.c_int
    lib.cuDeviceGetCount.argtypes = [ctypes.POINTER(ctypes.c_int)]
    lib.cuDeviceGetCount.restype = ctypes.c_int
    count = ctypes.c_int()
    if lib.cuInit(0) != 0:
        return 0
    if lib.cuDeviceGetCount(ctypes.byref(count)) != 0:
        return 0
    return count.value


def cuda_available() -> bool:
    """Return whether a CUDA device is available.

    This probe must be safe on machines without NVIDIA/CUDA because
    test decorators call it during pytest collection.
    """
    try:
        return cuda_device_count() > 0
    except (RuntimeError, OSError):
        return False
