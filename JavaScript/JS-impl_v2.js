
const boundCon = new Enumerator([periodic, fixed])

class Mesh {
    constructor(Nx, Ny, dx, dy) {
        this.Nx = Nx;
        this.Ny = Ny;
        this.dx = dx;
        this.dy = dy;
        this.inv_dxdy = 1 / (dx * dy);
        this.boundCons = {
            "xBoundL": boundCon.periodic,
            "xBoundR": boundCon.periodic,
            "yBoundD": boundCon.periodic,
            "yBoundU": boundCon.periodic
        }
        this.data = new Array((Nx + 2) * (Ny + 2)).fill(0);
        this.bounds = {
            "xBoundL": new Array(Ny).map((_, idx) => data[idx * Nx + 0.]),
            "xBoundR": new Array(Ny).map((_, idx) => data[idx * Nx + Nx]),
            "yBoundD": new Array(Ny).map((_, idx) => data[0. * Nx + idx]),
            "yBoundU": new Array(Ny).map((_, idx) => data[Ny * Nx + idx]),
        }
    }
    assign_boundary(xBoundL = boundCon.periodic, xBoundR = boundCon.periodic, yBoundD = boundCon.periodic, yBoundU = boundCon.periodic) {
        this.boundCons.xBoundL = xBoundL;
        this.boundCons.xBoundR = xBoundR;
        this.boundCons.yBoundD = yBoundD;
        this.boundCons.yBoundU = yBoundU;
    }
    get(x, y) {
        return this.data[(y + 1) * Nx + x + 1];
    }
    init_rand_value(baseValue, minVar, maxVar) {
        this.data.map(
            (_, _) => baseValue + (maxVar - minVar) * Math.random() + minVar
        );
    }
    update_mesh(ker_fun) {
        const new_mesh = new Array((this.Nx + 2) * (this.Ny + 2)).fill(0);
        for (let j = 1; j < this.Ny + 1; j++) {
            for (let i = 1; i < this.Nx + 1; i++) {
                const left = j * this.Nx + i - 1;
                const right = j * this.Nx + i + 1;
                const down = (j - 1) * this.Nx + i;
                const up = (j + 1) * this.Nx + i;
                const center = j * this.Nx + i;
                new_mesh[j * Nx + i] = ker_fun()
            }
        }
        for (let bnd in this.boundCons) {
            if (this.boundCons[bnd] === boundCon.periodic) {
                this.bounds[bnd].map()
            }
            if (bnd === boundCon.fixed) {

            }
        }
    }
    lap_ker() {
        return this.inv_dxdy * (
            this.data[left] +
            this.data[right] +
            this.data[up] +
            this.data[down] -
            4 * this.data[center]
        );
    }
}

class CH_Solver {
    constructor(M, dbulk_dc, dint_dc, others = undefined) {
        this.M = M
        this.dbulk_dc = dbulk_dc;
        this.dint_dc = dint_dc;
        if (others !== undefined) {
            this.others = others
        }
    }
    dF_dc() {
        return ((this.dbulk_dc + this.dint_dc).laplacian() * this.M).laplacian()
    }
}