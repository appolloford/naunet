from naunet.reactions.converter import ExpressionConverter


def test_c_converter():

    converter = ExpressionConverter("C")
    converter.read("pow(x, 2) * y[IDX_H2]")
    assert f"{converter}" == "pow(x, 2) * y[IDX_H2]"
    assert f"{converter:fortran}" == "x**2 * n(idx_H2)"

    converter.read("sqrt( exp(pow(x, 2)) )")
    assert f"{converter}" == "sqrt(exp(pow(x, 2)))"
    assert f"{converter:fortran}" == "sqrt(exp(x**2))"


def test_fortran_converter():

    converter = ExpressionConverter("Fortran")
    converter.read("x**2*sqrt(z)*n(idx_H2)")
    assert f"{converter}" == "x**2 * sqrt(z) * n(IDX_H2)"
    assert f"{converter:c}" == "pow(x, 2) * sqrt(z) * y[IDX_H2]"