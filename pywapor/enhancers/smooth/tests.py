import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from .whittaker import whittaker, cross_val_lmbda, wt1, wt2, second_order_diff_matrix, cve1, cve2

def cross_val_lmbda_ref(y, x, lmbdas = np.logspace(-2, 3, 100)):
    # Determine vector length.
    m = len(x)
    # Create arrays of m timeseries in which each time one value has been removed.
    Y_ = np.where(np.eye(m), np.nan, np.repeat(y[np.newaxis, :], m, axis = 0))
    # Smooth the timeeseries and make an estimate for the removed value for each lambda.
    z = whittaker(Y_, x, lmbdas = lmbdas)
    # Retrieve estimated timeseries for each lambda.
    y_hat = np.diagonal(z, axis1 = 0, axis2 = 2)
    # Calculate the cross validation standard error for each lambda.
    cves = np.sqrt(np.nanmean(((y - y_hat)**2), axis = 1))
    # Select the lambda for which the error is smallest.
    lmbda_sel = lmbdas[np.argmin(cves)]
    return lmbda_sel, cves

def test_cross_val(y, x):
    lmbdas = np.logspace(-2, 2, 20)
    lmbda_sel1, cves1 = cross_val_lmbda_ref(y, x, lmbdas = lmbdas)
    lmbda_sel2, cves2 = cross_val_lmbda(y, x, lmbdas = lmbdas)
    # assert np.isclose(lmbda_sel1, lmbda_sel2)
    plt.plot(lmbdas, cves1, label = "1")
    plt.plot(lmbdas, cves2, label = "2")
    plt.legend()

def test_shapes_z(y, x, xn = 2, yn = 5, n = 4, m = 3):
    t = len(y)

    Y1 = y[:]
    Y2 = np.repeat(y, xn*yn).reshape((xn,yn,t), order = "F")
    Y3 = np.repeat(y, n).reshape((n,t), order = "F")

    LMB1 = 100.
    LMB2 = np.logspace(1,3,m)
    LMB3 = np.logspace(1,3,xn*yn).reshape((xn, yn))

    D = second_order_diff_matrix(x, normalize = True)

    z11 = wt1(Y1, D, LMB1)
    assert z11.shape == (t,)

    z21 = wt1(Y2, D, LMB1)
    assert z21.shape == (xn,yn,t)

    z31 = wt1(Y3, D, LMB1)
    assert z31.shape == (n,t)

    z12 = wt1(Y1, D, LMB2) #
    assert z12.shape == (m,t)

    z22 = wt2(Y2, D, LMB2) #
    assert z22.shape == (xn,yn,m,t)

    z32 = wt2(Y3, D, LMB2) #
    assert z32.shape == (n,m,t)

    z13 = wt1(Y1, D, LMB3) #
    assert z13.shape == (xn,yn,t)

    z23 = wt1(Y2, D, LMB3) #!
    assert z23.shape == (xn,yn,t)

def test_shapes_cve(y, x, xn = 2, yn = 5, n = 4, m = 3):
    t = len(y)

    Y1 = y[:]
    Y2 = np.repeat(y, xn*yn).reshape((xn,yn,t), order = "F")
    Y3 = np.repeat(y, n).reshape((n,t), order = "F")

    LMB1 = 100.
    LMB2 = np.logspace(1,3,m)
    LMB3 = np.logspace(1,3,xn*yn).reshape((xn, yn))

    D = second_order_diff_matrix(x, normalize = True)

    cve11 = cve1(LMB1, Y1, D)
    assert cve11.shape == ()

    cve21 = cve1(LMB1, Y2, D)
    assert cve21.shape == (xn,yn)

    cve31 = cve1(LMB1, Y3, D)
    assert cve31.shape == (n,)

    cve12 = cve1(LMB2, Y1, D)
    assert cve12.shape == (m,)

    cve22 = cve2(LMB2, Y2, D)
    assert cve22.shape == (xn,yn,m)

    cve32 = cve2(LMB2, Y3, D)
    assert cve32.shape == (n,m)

    cve13 = cve1(LMB3, Y1, D)
    assert cve13.shape == (xn,yn)

    cve23 = cve1(LMB3, Y2, D)
    assert cve23.shape == (xn,yn)

def test_whittaker_main(ds, x, xn = 2, yn = 5, n = 4, m = 3):

    t = len(ds.time)

    base = 10

    Y1 = ds.isel(y = base, x = base)["ndvi"]
    Y2 = ds.isel(y = slice(base, base + yn), x = slice(base, base + xn))["ndvi"].transpose("y","x","time")
    Y3 = ds.isel(y = slice(base, base + n), x = base).rename({"y": "n"})["ndvi"]

    LMB1 = xr.DataArray(100.)
    LMB2 = xr.DataArray(np.logspace(1,3,m), dims = ["lmbda"], coords = {"lmbda": np.logspace(1,3,m)})
    LMB3 = xr.DataArray(np.logspace(1,3,xn*yn).reshape((xn, yn)), dims = ["x", "y"], coords = {"x": Y2.x, "y": Y2.y}).transpose("y","x")
    
    out1 = whittaker(Y1, x, lmbdas = LMB1, axis = -1)
    assert isinstance(out1, xr.DataArray)
    assert np.all([v == {"time": t}[k] for k, v in out1.sizes.items()])

    out2 = whittaker(Y2, x, lmbdas = LMB1, axis = -1)
    assert isinstance(out2, xr.DataArray)
    assert np.all([v == {"x": xn, "y": yn, "time": t}[k] for k, v in out2.sizes.items()])

    out3 = whittaker(Y3, x, lmbdas = LMB1, axis = -1)
    assert isinstance(out3, xr.DataArray)
    assert np.all([v == {"n": n, "time": t}[k] for k, v in out3.sizes.items()])

    out4 = whittaker(Y1, x, lmbdas = LMB2, axis = -1)
    assert isinstance(out4, xr.DataArray)
    assert np.all([v == {"lmbda": m, "time": t}[k] for k, v in out4.sizes.items()])

    out5 = whittaker(Y2, x, lmbdas = LMB2, axis = -1)
    assert isinstance(out5, xr.DataArray)
    assert np.all([v == {"lmbda": m, "time": t, "x": xn, "y": yn}[k] for k, v in out5.sizes.items()])

    out6 = whittaker(Y3, x, lmbdas = LMB2, axis = -1)
    assert isinstance(out6, xr.DataArray)
    assert np.all([v == {"lmbda": m, "time": t, "n": n}[k] for k, v in out6.sizes.items()])

    out7 = whittaker(Y1, x, lmbdas = LMB3, axis = -1)
    assert isinstance(out7, xr.DataArray)
    assert np.all([v == {"x": xn, "time": t, "y": yn}[k] for k, v in out7.sizes.items()])

    out8 = whittaker(Y2, x, lmbdas = LMB3, axis = -1)
    assert isinstance(out8, xr.DataArray)
    assert np.all([v == {"x": xn, "time": t, "y": yn}[k] for k, v in out8.sizes.items()])

    ######

    out1 = whittaker(Y1.values, x, lmbdas = LMB1.values, axis = -1)
    assert out1.shape == (t,)

    out2 = whittaker(Y2.values, x, lmbdas = LMB1.values, axis = -1)
    assert out2.shape == (yn,xn,t)

    out3 = whittaker(Y3.values, x, lmbdas = LMB1.values, axis = -1)
    assert out3.shape == (n,t)

    out4 = whittaker(Y1.values, x, lmbdas = LMB2.values, axis = -1)
    assert out4.shape == (m,t)

    out5 = whittaker(Y2.values, x, lmbdas = LMB2.values, axis = -1)
    assert out5.shape == (yn,xn,m,t)

    out6 = whittaker(Y3.values, x, lmbdas = LMB2.values, axis = -1)
    assert out6.shape == (n,m, t)

    out7 = whittaker(Y1.values, x, lmbdas = LMB3.values, axis = -1)
    assert out7.shape == (yn,xn,t)

    out8 = whittaker(Y2.values, x, lmbdas = LMB3.transpose("y","x").values, axis = -1)
    assert out8.shape == (yn,xn,t)

def test_cve_main(ds, x, xn = 2, yn = 5, n = 4, m = 3):

    t = len(ds.time)

    base = 10

    Y1 = ds.isel(y = base, x = base)["ndvi"]
    Y2 = ds.isel(y = slice(base, base + yn), x = slice(base, base + xn))["ndvi"].transpose("y","x","time")
    Y3 = ds.isel(y = slice(base, base + n), x = base).rename({"y": "n"})["ndvi"]

    LMB1 = xr.DataArray(100.)
    LMB2 = xr.DataArray(np.logspace(1,3,m), dims = ["lmbda"], coords = {"lmbda": np.logspace(1,3,m)})
    LMB3 = xr.DataArray(np.logspace(1,3,xn*yn).reshape((xn, yn)), dims = ["x", "y"], coords = {"x": Y2.x, "y": Y2.y}).transpose("y","x")
    
    out1, cves1 = cross_val_lmbda(Y1, x, lmbdas = LMB1)
    assert isinstance(out1, type(None))
    assert cves1.ndim == 0

    out2, cves2 = cross_val_lmbda(Y2, x, lmbdas = LMB1)
    assert isinstance(out2, type(None))
    assert np.all([v == {"x": xn, "y": yn,}[k] for k, v in cves2.sizes.items()])

    out3, cves3 = cross_val_lmbda(Y3, x, lmbdas = LMB1)
    assert isinstance(out3, type(None))
    assert np.all([v == {"n": n}[k] for k, v in cves3.sizes.items()])

    out4, cves4 = cross_val_lmbda(Y1, x, lmbdas = LMB2)
    assert np.all([v == {"lmbda": m}[k] for k, v in cves4.sizes.items()])
    assert out4.ndim == 0

    out5, cves5 = cross_val_lmbda(Y2, x, lmbdas = LMB2)
    assert np.all([v == {"lmbda": m, "x": xn, "y": yn}[k] for k, v in cves5.sizes.items()])
    assert np.all([v == {"x": xn, "y": yn}[k] for k, v in out5.sizes.items()])

    out6, cves6 = cross_val_lmbda(Y3, x, lmbdas = LMB2)
    assert np.all([v == {"lmbda": m, "n": n}[k] for k, v in cves6.sizes.items()])
    assert np.all([v == {"n": n}[k] for k, v in out6.sizes.items()])

    out7, cves7 = cross_val_lmbda(Y1, x, lmbdas = LMB3)
    assert isinstance(out7, type(None))   
    assert np.all([v == {"x": xn, "y": yn}[k] for k, v in cves7.sizes.items()])

    out8, cves8 = cross_val_lmbda(Y2, x, lmbdas = LMB3)
    assert isinstance(out8, type(None))   
    assert np.all([v == {"x": xn, "y": yn}[k] for k, v in cves8.sizes.items()])

    ######

    out1, cves1 = cross_val_lmbda(Y1.values, x, lmbdas = LMB1.values)
    assert isinstance(out1, type(None))

    out2, cves2 = cross_val_lmbda(Y2.values, x, lmbdas = LMB1.values)
    assert isinstance(out2, type(None))

    out3, cves3 = cross_val_lmbda(Y3.values, x, lmbdas = LMB1.values)
    assert isinstance(out3, type(None))

    out4, cves4 = cross_val_lmbda(Y1.values, x, lmbdas = LMB2.values)
    # assert isinstance(out3, type(None))   

    out5, cves5 = cross_val_lmbda(Y2.values, x, lmbdas = LMB2.values)
    # assert isinstance(out3, type(None))   

    out6, cves6 = cross_val_lmbda(Y3.values, x, lmbdas = LMB2.values)
    # assert isinstance(out3, type(None))   

    out7, cves7 = cross_val_lmbda(Y1.values, x, lmbdas = LMB3.values)
    assert isinstance(out7, type(None))   

    out8, cves8 = cross_val_lmbda(Y2.values, x, lmbdas = LMB3.values)
    assert isinstance(out8, type(None))  

    
