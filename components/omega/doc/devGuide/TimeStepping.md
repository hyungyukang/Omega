(omega-dev-time-stepping)=

# Time stepping

Omega supports choice of multiple time stepping schemes. The time stepping
module provides implementations of different time stepping schemes and general
infrastructure for managing time steppers. To achieve modularity, each time
stepper is implemented as a separate class. To ensure common interface and
enable sharing of common code and data, every time stepper inherits from the
`TimeStepper` base class.

## Time stepper type
An enumeration listing all implemented schemes is provided. It needs to be extended every time a new time stepper
is added. It is used to identify a time stepper at run time.
```c++
enum class TimeStepperType {
   ForwardBackward,
   RungeKutta4,
   RungeKutta2,
   SplitExplicitRK2,
   UnsplitRK2,
   Sundials,
   Invalid
};
```

## TimeStepper base class
The `TimeStepper` class is a base class for all Omega time steppers. It stores
as members common data needed by every time stepper. These include a name and
time stepping method as well as some time tracking variables, including
a Clock, TimeStep, StartTime, StopTime (optional) and an EndAlarm (optional).
The StopTime and EndAlarm are optional, because in coupled E3SM simulations the
components do not know the StopTime and therefore no EndAlarm can be created.
It also stores a StepCount, the number of `doStep` calls made on the
instance since creation.
Its `doStep` method defines the interface that every time stepper needs to
implement. In addition to this method, it provides a number of other methods
that can divided into two groups. The first group is for general time stepper
management. The second group is utility functions that provide common
functionality that can be re-used in implementation of different time steppers.

### Do step method
The `doStep` method is the main method of every time stepper class. It exists
in the base class as a pure virtual method, meaning it defines an interface but
has no implementation. Its signature is
```c++
    void TimeStepper::doStep(OceanState *State, TimeInstant &SimTime) const;
```
It is meant to advance `State` by one time step, from `Time` to `Time +
TimeStep`. Classes implementing concrete time stepping schemes need to provide
implementations of this method. As part of the time step, this routine advances
the Clock associated with the time stepper and returns the current simulation
time. For the default time stepper, this is advances the primary model clock.
Each implementation also increments StepCount, a member variable that
persists on the instance across calls, so it can be used to identify the
first `doStep` call even when `doStep` is invoked repeatedly across
separate calls to `ocnRun`, as is done by the coupled driver.

### Time stepper management methods

#### Initialization

Two static methods:
```c++
TimeStepper::init1();
TimeStepper::init2();
```
initialize the default time stepper in two phases. The first phase must be
called early in the initialization process. It creates and time stepper
instance, fills it with various time quantities (StartTime, EndTime, TimeStep)
and creates a Clock (the Model Clock for the simulation) and an EndAlarm to
ring when this portion of the simulation should stop. The second phase is
called after the default tendency, auxiliary state, horizontal mesh, and halo
objects have been initialized. It attaches pointers to these objects for later
use during the time stop phase. A pointer to the default time stepper can be
retrieved at any time using:
```c++
TimeStepper* DefTimeStepper = TimeStepper::getDefault();
```

In standalone simulations, `TimeStepper::init1();` will read the StartTime from
the configuration file. For coupled simulations the StartTime is provided by
the coupler, overriding the value defined in the configuration file.
In this case, the first phase of time stepper initialization would look like:
```c++
// initialize TimeStepper with coupler provided start type and start time
TimeStepper::init1(StartType, StartTime);
```

#### Creation of non-default time steppers

A non-default time stepper can be created from a string `Name`, time stepper
type `Type`, `TimeStep`, `StartTime`, `EndTime`, tendencies `Tend`,
auxiliary state `AuxState`, horizontal mesh `Mesh`, and halo layer `MyHalo`
```c++
TimeStepper*  NewTimeStepper = TimeStepper::create(Name, Type,
        Tend, AuxState, Mesh, MyHalo, TimeStep, StartType, StartTime,
        StopType, EndTime, Duration);
```
where the EndTime and Duration are optional arguments depending on StopType.

The TimeStepper constructors are filled initially with a StopTime based on
the simulation CurrentTime. In many cases (including the default time stepper),
the constructor is called before the CurrentTime is reset on restart. For the
StartType Continue option, the StopTime and EndAlarm must be reset by calling
the ``MyStepper->resetEndAlarm()`` function once the current simulation time
has been updated from the restart metadata.

For convenience, the create routine returns a pointer to the newly created
time stepper. Given its name, a pointer to a named time stepper can be obtained
at any time by calling the static `get` method:
```c++
TimeStepper* NewTimeStepper = TimeStepper::get(Name);
```

#### Changing the time step
The time step is set on creation, but if the time step must change, the
function
```c++
Stepper->changeTimeStep(TimeStep);
```
can be used where `TimeStep` is an instance of `TimeInterval` class.

#### Getters
Given a pointer to a `TimeStepper` you can obtain various class members.
Pointers to the associated Clock and EndAlarm can also be retrieved.
```c++
TimeStepperType Type = Stepper->getType();
int NTimeLevels = Stepper->getNTimeLevels();
std::string Name = Stepper->getName();
TimeInterval TimeStep = Stepper->getTimeStep();
TimeStepperStartType StartType = Stepper->getStartType();
TimeInstant StartTime = Stepper->getStartTime();
TimeStepperStopType StopType = Stepper->getStopType();
TimeInstant StopTime = Stepper->getStopTime();
TimeInterval Duration = Stepper->getDuration();
Clock *ModelClock = Stepper->getClock();
Alarm *EndAlarm = Stepper->getEndAlarm();
```
Note that StopTime, Duration and EndAlarm are set to appropriate values
no matter what the StopType is (though for the Coupled option, they are
set far into the future).

#### Removal of time steppers
To erase a specific named time stepper use `erase`
```c++
TimeStepper::erase(Name);
```
To clear all instances do:
```c++
TimeStepper::clear();
```

### Utility methods

In addition to the above management methods, the base class provides common
functions that are useful in implementing different time steppers. Given two
pointers to `OceanState` objects `State1` and `State2`, two integers
representing two time levels `TimeLevel1` and `TimeLevel2`, and coefficient
`Coeff` of type `TimeInterval` calling
```c++
updateStateByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
conceptually performs the following update
`State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend`
element-wise for every state variable. To perform this update only for
pseudo-thickness do
```c++
updateThicknessByTend(State1, TimeLevel1, *State2, TimeLevel2, Coeff);
```
and only for normal velocity do
```c++
updateVelocityByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
Similarly for tracers, the method for updating tracers is
```c++
updateTracersByTend(NextTracers, CurTracers, State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
Additionally, for the RK4 method, there are methods for initializing, accumulating, and
finalizing the tracer integration for a given time step
```c++
weightTracers(NextTracers, CurTracers, CurState, TimeLevel1);
accumulateTracersUpdate(AccumTracer, Coeff);
finalizeTracersUpdate(NextTracers, State, TimeLevel);
```

## Implemented time steppers
The following time steppers are currently implemented
| Class name | Enum value | Scheme |
| ---------- | -------------------------- | ------ |
| ForwardBackwardStepper | TimeStepperType::ForwardBackward | forward-backward |
| RungeKutta2Stepper | TimeStepperType::RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4Stepper | TimeStepperType::RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
| SundialsStepper | TimeStepperType::Sundials | explicit Runge-Kutta method chosen in the configuration, taken by SUNDIALS ARKODE |

### SundialsStepper
The `SundialsStepper` takes each step with the ERKStep integrator of the
SUNDIALS ARKODE package, so any explicit Butcher table can be used. The table
is loaded when the stepper is constructed, and the stage fractions `c` and the
projection fractions passed to `computeAllTendencies` as `ProjDt` (the next
stage's `c`, and 1 for the last stage, as in `RungeKutta4Stepper`) are derived
from it. The ARKODE integrator is created on the first call to `doStep`.

The state ARKODE advances is an MPI ManyVector of three pieces: pseudo-thickness,
normal velocity and thickness-weighted tracers. Each piece is an Omega array
wrapped, without copying, as a SUNDIALS Kokkos vector (`SundialsVector.h`), so
SUNDIALS kernels run in Omega's execution space on every backend. The wrappers
are created once per array, because Omega rotates time levels without moving
arrays. `SundialsVector.h` also gives the Kokkos vector a fused linear
combination, which forms each stage in one kernel and adds the terms in the
order of Omega's own Runge-Kutta updates.

Each `doStep`
1. weights the current tracers by thickness into the next tracer time level,
   prescribes the state and computes the KPP fields, as `RungeKutta4Stepper`
   does before its first stage
2. resets ARKODE to the current time level at `t = 0` and takes one fixed step
   of length `TimeStep` into a vector wrapping the next time level, so ARKODE
   forms every intermediate stage directly in Omega's arrays. The first stage
   is evaluated on the current time level. Before each later stage, a stage
   postprocessing function masks the velocity, prescribes the state at the
   stage time, converts the tracers to concentrations and exchanges halos. The
   right-hand side function computes the tendencies and copies them to ARKODE's
   vector, with zeros outside the active layers so ARKODE leaves those values
   unchanged.
3. converts the tracers back to concentrations, rotates the time levels,
   applies implicit vertical mixing, validates the state and advances the clock

ARKODE allocates its work vectors (one per stage plus about six more) on the
first steps only. In the time stepper unit test, `RK4` and `RK2` reproduce
`RungeKutta4Stepper` and `RungeKutta2Stepper` to round-off; the test checks
both, and the convergence order of several methods. With a prescribed thickness
or velocity, `RK4` differs from `RungeKutta4Stepper`, which prescribes its
later stages at the accumulated times `t + sum(RKA) dt` rather than at
`t + c dt`.

SUNDIALS is built in the precision of Omega's `Real`. The single-precision copy
of the Omega library built for some tests cannot share vectors with it, so
there the constructor aborts.
