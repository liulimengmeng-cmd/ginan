#include "../src/cpp/pea/pppArJointGate.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

int main(int argc, char** argv)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    if (argc == 2)
    {
        std::ifstream input(argv[1]);
        std::string line;
        std::getline(input, line);
        require(line == "GINAN_AR_SNAPSHOT_V1", "snapshot schema");
        std::getline(input, line);
        MatrixXd matrices[6];
        for (int i = 0; i < 6; ++i)
        {
            std::string name;
            int rows, cols;
            input >> name >> rows >> cols;
            matrices[i].resize(rows, cols);
            for (int row = 0; row < rows; ++row)
                for (int col = 0; col < cols; ++col)
                    input >> matrices[i](row, col);
        }
        require(bool(input), "snapshot matrix parse");
        const auto S = (matrices[2] * matrices[1] * matrices[2].transpose()
                        + matrices[4]).eval();
        auto result = pppArJointGate(S, matrices[3].col(0), 4);
        std::cout << std::setprecision(17) << result.status << " " << result.nis
                  << " " << result.threshold << " " << result.alpha << "\n";
        return 0;
    }
    MatrixXd S(2, 2);
    S << 1, .999, .999, 1;
    VectorXd v(2);
    v << .2, -.2;
    const MatrixXd before = S;
    const VectorXd vBefore = v;
    require(!pppArJointGate(S, v, 4).passed, "joint conflict must reject");
    require(pppArJointGate(S.topLeftCorner(1, 1), v.head(1), 4).passed,
            "marginal candidate must pass synthetic example");
    v << .2, .2;
    require(pppArJointGate(S, v, 4).passed, "compatible common error must pass");
    v = vBefore;
    require(!pppArJointGate(S * 100, v * 10, 4).passed, "unit invariance");
    require(S == before && v == vBefore, "gate must not mutate inputs");
    S(0, 1) = .5;
    require(std::string(pppArJointGate(S, v, 4).status) == "ASYMMETRIC_COVARIANCE",
            "asymmetry must fail closed");
    S.setOnes();
    require(!pppArJointGate(S, v, 4).passed, "singular covariance must reject");
    S.setIdentity();
    v[0] = std::numeric_limits<double>::quiet_NaN();
    require(!pppArJointGate(S, v, 4).passed, "nonfinite innovation must reject");
    v.setZero();
    require(!pppArJointGate(S, v, 0).passed, "invalid sigma must reject");
    require(!pppArJointGate(MatrixXd(), VectorXd(), 4).passed, "empty block must reject");
    MatrixXd one = MatrixXd::Identity(1, 1);
    VectorXd near(1);
    near << 3.999;
    auto pass = pppArJointGate(one, near, 4);
    require(pass.passed && std::abs(pass.threshold - 16) < 1e-9,
            "one dimensional four sigma reference");
    near << 4.001;
    require(!pppArJointGate(one, near, 4).passed, "outside reference must reject");
    std::cout << "PASS: correlated rejection, compatible pass, units, immutability, "
              << "invalid inputs, and analytic one-dimensional boundary\n";
}
