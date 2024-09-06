from naunet.reactions.converter import ExpressionConverter


def test_c_converter():
    converter = ExpressionConverter("C")
    converter.read("pow(x, 2) * y[IDX_H2]")
    assert f"{converter}" == "pow(x, 2) * y[IDX_H2]"
    assert f"{converter:fortran}" == "(x)**2 * n(idx_H2)"

    converter.read("sqrt( exp(pow(x, 2)) )")
    assert f"{converter}" == "sqrt(exp(pow(x, 2)))"
    assert f"{converter:fortran}" == "sqrt(exp((x)**2))"

    converter.read("8.16e-10 * pow(Tgas/300.0, 0.0) * exp(-164.9/Tgas)")
    assert f"{converter}" == "8.16e-10 * pow(Tgas/300.0, 0.0) * exp(-164.9/Tgas)"
    assert f"{converter:fortran}" == "8.16e-10 * (Tgas/300.0)**0.0 * exp(-164.9/Tgas)"

    converter.read("8.16e-10 * pow(xr + cr, 2.0) * exp(-164.9/Tgas)")
    assert f"{converter}" == "8.16e-10 * pow(xr + cr, 2.0) * exp(-164.9/Tgas)"
    assert f"{converter:fortran}" == "8.16e-10 * (xr + cr)**2.0 * exp(-164.9/Tgas)"


def test_fortran_converter():
    converter = ExpressionConverter("Fortran")
    converter.read("x**2*sqrt(z)*n(idx_H2)")
    assert f"{converter}" == "x**2 * sqrt(z) * n(IDX_H2)"
    assert f"{converter:c}" == "pow(x, 2) * sqrt(z) * y[IDX_H2]"
